"""
Author: Junkai Zhang (Kain)
Created:  '2025-09-15'
Last Modified:2025-09-15
Email: junkaiz@vt.edu
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import LogFormatterSciNotation, LogLocator
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import FuncFormatter


def load_hic_matrix(file_path, delimiter=','):
    """
    Load a Hi-C file and return a 2D matrix.
    If the file delimiter is not a comma, change the delimiter parameter.
    """
    return np.loadtxt(file_path, delimiter=delimiter)

def compute_row_normalized_P(matrix, exclude_offset):
    """
    Compute the row-normalized P matrix:
      1. Define mask_{ij} = (|i-j| > exclude_offset)
      2. A = matrix.copy(); A[~mask] = 0
      3. row_sum[i] = Σ_k A[i,k]; P[i,j] = A[i,j] / row_sum[i]
    """
    N = matrix.shape[0]
    I, J = np.indices((N, N))
    mask = (np.abs(I - J) > exclude_offset)
    A = matrix.copy()
    A[~mask] = 0
    row_sums = A.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return A / row_sums


def plot_difference(P_wt, P_lm, exclude_offset):
    """
    Plot the difference map D = P_lm - P_wt.
    Set the main diagonal and ±exclude_offset to 0,
    and fix the red/blue cutoff at ±1e-3.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # 1. Compute difference and mask diagonal band
    D = P_lm - P_wt
    N = D.shape[0]
    I, J = np.indices((N, N))
    D[np.abs(I - J) <= exclude_offset] = 0

    # 2. Fix cutoff at ±1e-3 to highlight small differences
    cutoff = 1e-3
    norm = TwoSlopeNorm(vcenter=0, vmin=-cutoff, vmax=cutoff)

    # 3. Mask zero values and plot
    D_masked = np.ma.masked_where(D == 0, D)
    fig, ax = plt.subplots(figsize=(6, 6), facecolor='white')
    img = ax.imshow(
        D_masked,
        cmap='RdBu_r',
        norm=norm,
        interpolation='nearest'
    )

    # 4. Beautify
    ax.set_title(f"Difference WT3h – WT_Rejuv (K={exclude_offset})", fontsize=12, pad=10)
    ax.set_xlabel('j'); ax.set_ylabel('i')
    ax.set_xticks([0, N//2, N-1]); ax.set_yticks([0, N//2, N-1])
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3)

    # 5. Configure colorbar with scientific notation
    cbar = fig.colorbar(
        img,
        ax=ax,
        shrink=0.8,
        format=FuncFormatter(lambda x, pos: f"{x:.1e}")
    )
    cbar.set_label("ΔP = P_LM – P_WT", rotation=270, labelpad=15)

    plt.tight_layout()
    plt.savefig("WT3h_WT(rejuv)_diff_enhanced.png", bbox_inches='tight')
    plt.show()


def main():
    # # === User input ===
    # wt_path = "Week14\hicmap\WT3h\WT3h.map"  # WT Hi-C matrix file path
    # lm_path = "Week14\hicmap\WT_after_dead3h\WT_after_dead3h.map"   # LM Hi-C matrix file path
    # delimiter = ','                      # Modify according to file delimiter
    # exclude_offset = 4                   # Exclude main diagonal and ±4
    # # ==================

    # # 1. Load and normalize
    # mat_wt = load_hic_matrix(wt_path, delimiter)
    # mat_lm = load_hic_matrix(lm_path, delimiter)
    # P_wt = compute_row_normalized_P(mat_wt, exclude_offset)
    # P_lm = compute_row_normalized_P(mat_lm, exclude_offset)

    # N = P_wt.shape[0]
    # # 2. Construct combined matrix: lower triangle WT, upper triangle LM
    # P_combo = np.zeros_like(P_wt)
    # tril_i, tril_j = np.tril_indices(N, k=-1)
    # triu_i, triu_j = np.triu_indices(N, k=1)
    # P_combo[tril_i, tril_j] = P_wt[tril_i, tril_j]
    # P_combo[triu_i, triu_j] = P_lm[triu_i, triu_j]
    # # Mask main diagonal and ±exclude_offset
    # I, J = np.indices((N, N))
    # mask_diag = (np.abs(I - J) <= exclude_offset)
    # P_combo[mask_diag] = 0

    # # 3. Plot
    # P_masked = np.ma.masked_where(P_combo == 0, P_combo)
    # cmap = LinearSegmentedColormap.from_list('custom_red', ['white', 'red'])
    # cmap.set_bad('white')

    # fig, ax = plt.subplots(figsize=(6, 6), facecolor='white')
    # nonzero = P_masked.compressed()
    # img = ax.imshow(
    #     P_masked,
    #     cmap=cmap,
    #     interpolation='nearest',
    #     norm=LogNorm(vmin=nonzero.min(), vmax=nonzero.max())
    # )
    # ax.set_title(f"WT3h (lower) & WT_Rejuv (upper) Hi-C P$_{{ij}}$ (log scale, K={exclude_offset})", fontsize=12, pad=10)
    # ax.set_xlabel('j')
    # ax.set_ylabel('i')
    # ax.set_xticks([0, N//2, N-1])
    # ax.set_yticks([0, N//2, N-1])
    # ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3)

    # cbar = fig.colorbar(
    #     img, ax=ax,
    #     label="P$_{ij}$ (log scale)",
    #     shrink=0.8,
    #     format=LogFormatterSciNotation(labelOnlyBase=True)
    # )
    # cbar.locator = LogLocator(base=10)
    # cbar.update_ticks()

    # plt.tight_layout()
    # plt.savefig("WT_vs_WT_triangle.png", bbox_inches='tight')
    # plt.show()
    #**************************************************************************************************
    #********************************* Below is the difference plot; above is the combined-triangle matrix plot *********************************
    #**************************************************************************************************
    # === User input ===
    wt_path = "Week14\hicmap\WT3h\WT3h.map"
    lm_path = "Week14\hicmap\WT_after_dead3h\WT_after_dead3h.map"
    delimiter = ','        # or '\t'
    exclude_offset = 4
    # ==================

    mat_wt = load_hic_matrix(wt_path, delimiter)
    mat_lm = load_hic_matrix(lm_path, delimiter)

    P_wt = compute_row_normalized_P(mat_wt, exclude_offset)
    P_lm = compute_row_normalized_P(mat_lm, exclude_offset)

    # Print entropy values (optional)
    # S_wt, spt_wt = compute_entropies(P_wt)
    # S_lm, spt_lm = compute_entropies(P_lm)
    # print(f"WT: S={S_wt:.4f}, S/N={spt_wt:.6f}")
    # print(f"LM: S={S_lm:.4f}, S/N={spt_lm:.6f}")
    D = P_lm - P_wt
    print("Min difference:", D.min())
    print("Max difference:", D.max())
    print("Mean difference:", D.mean())

    # Plot difference
    plot_difference(P_wt, P_lm, exclude_offset)

if __name__ == '__main__':
    main()
