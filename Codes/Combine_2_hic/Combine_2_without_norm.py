"""
Author: Junkai Zhang (Kain)
Created:  '2025-09-18'
Last Modified:2025-09-18
Email: junkaiz@vt.edu
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
from matplotlib.ticker import LogFormatterSciNotation, FuncFormatter

def load_hic_matrix(path, delimiter=','):
    """Load the raw Hi-C matrix without normalization."""
    return np.loadtxt(path, delimiter=delimiter)

def plot_triangle_raw(mat_wt, mat_lm, exclude_offset, out_png="WT_vs_LM_triangle_raw.png"):
    """
    In one plot:
      - Lower triangle uses WT raw matrix
      - Upper triangle uses LM raw matrix
      - Main diagonal and ±exclude_offset are masked as 0 (white)
      - Logarithmic colormap
    """
    N = mat_wt.shape[0]
    # Construct combined matrix
    M = np.zeros_like(mat_wt)
    tri_i, tri_j = np.tril_indices(N, k=-1)
    upr_i, upr_j = np.triu_indices(N, k=1)
    M[tri_i, tri_j] = mat_wt[tri_i, tri_j]
    M[upr_i, upr_j] = mat_lm[upr_i, upr_j]
    # Mask diagonal and ±K
    I, J = np.indices((N, N))
    M[np.abs(I - J) <= exclude_offset] = 0

    # Mask zero values as white
    M_masked = np.ma.masked_where(M == 0, M)
    cmap = LinearSegmentedColormap.from_list('custom_red', ['white','red'])
    cmap.set_bad('white')

    # Plot
    fig, ax = plt.subplots(figsize=(6,6), facecolor='white')
    # Use logarithmic scale (vmin/vmax taken from non-zero values automatically)
    vals = M_masked.compressed()
    img = ax.imshow(
        M_masked,
        cmap=cmap,
        norm=LogNorm(vmin=vals.min(), vmax=vals.max()),
        interpolation='nearest'
    )
    ax.set_title(f"WT (lower) & LM (upper) raw counts (K={exclude_offset})", pad=10)
    ax.set_xlabel('j'); ax.set_ylabel('i')
    ax.set_xticks([0, N//2, N-1]); ax.set_yticks([0, N//2, N-1])
    ax.grid(True, linestyle='--', alpha=0.3)

    # Colorbar
    cbar = fig.colorbar(
        img, ax=ax, shrink=0.8,
        format=LogFormatterSciNotation(labelOnlyBase=True)
    )
    cbar.set_label("Contact count (log scale)", rotation=270, labelpad=15)
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches='tight')
    plt.show()

def plot_difference_raw(mat_wt, mat_lm, exclude_offset, out_png="LM_minus_WT_diff_raw.png"):
    """
    Difference map D = mat_lm - mat_wt (raw counts):
      - Use TwoSlopeNorm to center at 0
      - Mask main diagonal and ±exclude_offset
      - Upper/lower limits set by 98th percentile of abs(D), 
        to emphasize common differences
    """

    # 1. Compute difference and mask diagonal band
    D = mat_lm - mat_wt
    N = D.shape[0]
    I, J = np.indices((N, N))
    D[np.abs(I - J) <= exclude_offset] = 0

    # 2. Mask zero values, set limits by 98th percentile of |D|
    D_masked = np.ma.masked_where(D == 0, D)
    abs_vals = np.abs(D_masked.compressed())
    clip = np.percentile(abs_vals, 98)    # 98th percentile
    norm = TwoSlopeNorm(vcenter=0, vmin=-clip, vmax=clip)

    # 3. Plot
    fig, ax = plt.subplots(figsize=(6, 6), facecolor='white')
    img = ax.imshow(
        D_masked,
        cmap='RdBu_r',
        norm=norm,
        interpolation='nearest'
    )
    ax.set_title(f"Difference LM – WT raw (K={exclude_offset})", fontsize=12, pad=10)
    ax.set_xlabel('j'); ax.set_ylabel('i')
    ax.set_xticks([0, N//2, N-1]); ax.set_yticks([0, N//2, N-1])
    ax.grid(True, linestyle='--', alpha=0.3)

    # 4. Colorbar
    cbar = fig.colorbar(
        img, ax=ax, shrink=0.8,
        format=FuncFormatter(lambda x, pos: f"{x:.1f}")
    )
    cbar.set_label("Δcount = LM − WT", rotation=270, labelpad=15)

    plt.tight_layout()
    plt.savefig(out_png, bbox_inches='tight')
    plt.show()

def main():
    # === User input ===
    wt_path = "Week14\hicmap\WT3h\WT3h.map"
    lm_path = "Week14\hicmap\WT_after_dead3h\WT_after_dead3h.map"
    delimiter     = ','      # or '\t'
    exclude_offset= 0
    # ==================

    mat_wt = load_hic_matrix(wt_path, delimiter)
    mat_lm = load_hic_matrix(lm_path, delimiter)
    D = mat_lm - mat_wt
    print("Min difference:", D.min())
    print("Max difference:", D.max())
    print("Mean difference:", D.mean())
    # Plot combined triangle (raw)
    plot_triangle_raw(mat_wt, mat_lm, exclude_offset)

    # Plot difference map (raw)
    plot_difference_raw(mat_wt, mat_lm, exclude_offset)

if __name__ == '__main__':
    main()
