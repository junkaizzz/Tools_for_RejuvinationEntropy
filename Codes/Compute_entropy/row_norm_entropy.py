"""
Author: Junkai Zhang (Kain)
Created:  '2025-09-15'
Last Modified:2025-09-15
Email: junkaiz@vt.edu
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import LogFormatterSciNotation, LogLocator
from pathlib import Path


def generate_test_matrix(N, regions, P1, P2):
    """
    Generate an N×N test Hi-C matrix:
      - Initialize all elements with P2 (inter-chromosomal contact probability)
      - For each chromosome interval, assign the submatrix within that region to P1 (intra-chromosomal contact probability)
    """
    matrix = np.full((N, N), P2)
    for (start, end) in regions:
        matrix[start:end, start:end] = P1
    return matrix


def create_chromosome_labels(N, regions, labels):
    """
    Generate label list based on chromosome intervals
    """
    chrom_labels = [''] * N
    for (region, label) in zip(regions, labels):
        start, end = region
        for i in range(start, end):
            chrom_labels[i] = label
    return chrom_labels


def compute_row_normalized_P(matrix, chrom_labels, exclude_offset):
    """
    1. Construct mask, excluding same chromosome pairs with |i-j| <= exclude_offset
    2. Set excluded values to zero and perform row normalization
    3. Return the normalized P matrix
    """
    N = matrix.shape[0]
    I, J = np.indices((N, N))
    chrom_arr = np.array(chrom_labels)
    same_chrom = (chrom_arr[:, None] == chrom_arr[None, :])
    intra_mask = np.abs(I - J) > exclude_offset
    mask = np.where(same_chrom, intra_mask, True)

    A = matrix.copy()
    A[~mask] = 0
    row_sums = A.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    P = A / row_sums
    return P


def compute_entropies_global(matrix, chrom_labels, exclude_offset):
    """
    Compute row-normalized entropies and return the P matrix as well
    """
    P = compute_row_normalized_P(matrix, chrom_labels, exclude_offset)
    logP = np.where(P > 0, np.log2(P), 0.0)
    S_i = -np.sum(P * logP, axis=1)
    total_entropy = float(np.sum(S_i))
    entropy_per_TAD = total_entropy / P.shape[0]
    return total_entropy, entropy_per_TAD, P


def main():
    # Parameters (modify directly)
    exclude_offset = 4  # Exclude main diagonal and ±1

    # Load real Hi-C matrix (please modify with actual file path and delimiter)
    file_path = "D:\Chro\HiC_0.2_evo10_av1min_av18.map"
    matrix = np.loadtxt(file_path, delimiter=',')

    # Build labels
    N = matrix.shape[0]
    regions = [(0, 225), (225, 438), (438, 659), (659, 966), (966, 985), (985, N)]
    labels = ["2L", "2R", "3L", "3R", "4", "X"]
    chrom_labels = create_chromosome_labels(N, regions, labels)

    # Compute entropy and P matrix
    total_entropy, entropy_per_TAD, P = compute_entropies_global(matrix, chrom_labels, exclude_offset)
    print(
    f"File: {file_path}\n"
    f"Exclude = {exclude_offset}\n"
    f"Total entropy (bits) = {total_entropy:.4f}\n"
    f"Entropy per TAD (bits) = {entropy_per_TAD:.6f}"
)

    # Plot: row-normalized P matrix, masked region (P==0) shown as white
    P_masked = np.ma.masked_where(P == 0, P)
    cmap = LinearSegmentedColormap.from_list('custom_red', ['white', 'red'])
    cmap.set_bad('white')
    fig, ax = plt.subplots(figsize=(6, 6), facecolor='white')
    nonzero = P_masked.compressed()
    img = ax.imshow(
        P_masked,
        cmap=cmap,
        interpolation='nearest',
        norm=LogNorm(vmin=0.0001, vmax=0.1)
    )
    ax.set_title(
        f"WT3h Hi-C matrix P$_{{ij}}$ (log scale, s$>{exclude_offset}$)",
        pad=10, fontsize=12
    )
    ax.set_xlabel('j')
    ax.set_ylabel('i')
    ax.set_xticks([0, N//2, N-1])
    ax.set_yticks([0, N//2, N-1])
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3)

    # Configure colorbar, show only scientific notation ticks, do not cover masked areas
    cbar = fig.colorbar(
        img, ax=ax,
        label="P$_{ij}$ (log scale)",
        shrink=0.8,
        format=LogFormatterSciNotation(labelOnlyBase=True)
    )
    cbar.locator = LogLocator(base=10)
    cbar.update_ticks()

    plt.tight_layout()
    fig.savefig("hic_matrix_heatmap.png", bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
