"""
Author: Junkai Zhang (Kain)
Created: 2025-09-18
Last Modified: 2025-09-18
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
    Generate an N×N synthetic Hi-C matrix:
      - Initialize all elements with P2 (inter-chromosomal contact probability)
      - For each chromosome interval, set the corresponding submatrix to P1 (intra-chromosomal contact probability)
    """
    matrix = np.full((N, N), P2)
    for (start, end) in regions:
        matrix[start:end, start:end] = P1
    return matrix


def create_chromosome_labels(N, regions, labels):
    """
    Generate chromosome labels for each bin index based on the given intervals.
    """
    chrom_labels = [''] * N
    for (region, label) in zip(regions, labels):
        start, end = region
        for i in range(start, end):
            chrom_labels[i] = label
    return chrom_labels


def compute_row_normalized_P(matrix, chrom_labels, exclude_offset):
    """
    1. Build a mask to exclude same-chromosome interactions when |i-j| <= exclude_offset.
    2. Zero out excluded entries and perform row-wise normalization.
    3. Return the normalized probability matrix P.
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
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    P = A / row_sums
    return P


def compute_entropies_global(matrix, chrom_labels, exclude_offset):
    """
    Compute row-normalized entropy for the whole matrix, and return both the entropy
    and the normalized probability matrix P.
    """
    P = compute_row_normalized_P(matrix, chrom_labels, exclude_offset)
    logP = np.where(P > 0, np.log2(P), 0.0)
    S_i = -np.sum(P * logP, axis=1)
    total_entropy = float(np.sum(S_i))
    entropy_per_TAD = total_entropy / P.shape[0]
    return total_entropy, entropy_per_TAD, P


def main():
    # Parameter settings (modify directly here)
    mode = 'block'      # 'uniform' or 'block'
    P1 = 1              # Intra-chromosomal probability
    P2 = 1e-2           # Inter-chromosomal probability
    exclude_offset = 1  # Exclude diagonal and ±1 bins

    # Build matrix and labels
    N = 1169
    regions = [(0,225),(225,438),(438,659),(659,966),(966,985),(985,N)]
    labels = ["2L","2R","3L","3R","4","X"]
    chrom_labels = create_chromosome_labels(N, regions, labels)

    # Generate synthetic matrix
    if mode == 'uniform':
        matrix = generate_test_matrix(N, regions, P1, P1)
    else:
        matrix = generate_test_matrix(N, regions, P1, P2)

    # Compute entropy and P matrix
    total_entropy, entropy_per_TAD, P = compute_entropies_global(matrix, chrom_labels, exclude_offset)
    print(f"Mode: {mode}\nP1 = {P1}, P2 = {P2}\nExclude = {exclude_offset}\n"
          f"Total entropy (bits) = {total_entropy:.4f}\n"
          f"Entropy per TAD (bits) = {entropy_per_TAD:.6f}\n")

    # Plot row-normalized P
    P_masked = np.ma.masked_where(P == 0, P)

    cmap = LinearSegmentedColormap.from_list('custom_red', ['white', 'red'])
    cmap.set_bad('white')

    fig, ax = plt.subplots(figsize=(6,6), facecolor='white')

    nonzero = P_masked.compressed()
    vmin, vmax = 0.0001, 0.1

    img = ax.imshow(
        P_masked,
        cmap=cmap,
        interpolation='nearest',
        norm=LogNorm(vmin=vmin, vmax=vmax)
    )

    ax.set_title(
        f"Row-normalized EPC HiC matrix P$_{{ij}}$ (log scale, K$={exclude_offset}$)",
        pad=10,
        fontsize=12
    )
    ax.set_xlabel('j')
    ax.set_ylabel('i')
    ax.set_xticks([0, 500, 1000])
    ax.set_yticks([0, 500, 1000])
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3)

    # Configure colorbar with scientific notation
    cbar = fig.colorbar(
        img, ax=ax,
        label="P$_{ij}$ (log scale)",
        shrink=0.8,
        format=LogFormatterSciNotation(labelOnlyBase=True)
    )
    cbar.locator = LogLocator(base=10)
    cbar.update_ticks()

    plt.tight_layout()
    fig.savefig(f"synthetic_{mode}_heatmap.png", bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
