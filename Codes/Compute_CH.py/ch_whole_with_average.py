# -*- coding: utf-8 -*-
"""
CH(s) from Hi-C probabilities using row-normalized P.
Outputs ONLY:
  - CH_comparison_WHOLE.png (whole-genome, across states)
  - CH_comparison_chr_<LABEL>.png for LABEL in [2L,2R,3L,3R,4,X] (across states)

Rules:
  - Teacher-approved row normalization: mask same-chrom |i-j|<=EXCLUDE_OFFSET; keep inter-chromosomal values
  - Use BOTH upper & lower triangles when forming P(s)
  - Y axis for WHOLE: 0–3 with label "CH(s) (×10⁻⁴)"
  - Running average smoothing:
        w = 100, centered:
        * use only full windows (no padding)
        * x-coordinate is the center of each window
        * plotted s-range becomes shorter after smoothing
  - X axis: genomic distance (not "bins"), same notation as Samira's graphs.

Author: Junkai(Kain) Zhang
Date: 2025-10-30
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ========= Basic Settings =========
DATA_ROOT = r"D:\Chro"                 # ← Change to your data root directory
OUT_DIR   = r"./out_ch_running_average_100final"
os.makedirs(OUT_DIR, exist_ok=True)

STATE_DIRS = {
    "WT_orig":        "WT3h_orig",
    "CD3h":           "CD3h",
    "LM3h":           "LM3h",
    "WT_after_CD3h":  "WT3h_after_CD3h",
    "WT_after_LM3h":  "WT3h_after_LM3h",
}

# Display names for the graph legend (use the same naming convention as entropy plots)
STATE_DISPLAY = {
    "WT_orig":        "WT(Original)",
    "CD3h":           "'Dead'",
    "LM3h":           "Lamins depleted",
    "WT_after_CD3h":  "Rejuvenated after 'Dead'",
    "WT_after_LM3h":  "Rejuvenated after Lamins depleted",
}

GLOB_PATTERN = "*.map"
DELIM = ','

EXCLUDE_OFFSET = 4
N_BINS  = 1169
REGIONS = [(0, 225), (225, 438), (438, 659), (659, 966), (966, 985), (985, N_BINS)]
LABELS  = ["2L", "2R", "3L", "3R", "4", "X"]

# ========= Genomic Distance Conversion =========
# Set the bin size consistent with Samira's notation (e.g., 10 kb per bin)
BIN_SIZE_BP = 10_000
X_UNIT = "Mb"   # Use "kb" if desired

def bins_to_distance(s_array):
    """Convert genomic distance (s in bins) into physical distance (kb/Mb)."""
    s_array = np.asarray(s_array, dtype=float)
    if X_UNIT.lower() == "mb":
        return s_array * BIN_SIZE_BP / 1e6
    elif X_UNIT.lower() == "kb":
        return s_array * BIN_SIZE_BP / 1e3
    else:
        return s_array * BIN_SIZE_BP / 1e6

# ========= Plot Aesthetics =========
plt.rcParams["figure.facecolor"] = "white"
COMBINED_FIGSIZE = (7.5, 5)

CH_SCALE = 1e-4
YMAX_WHOLE_UNITS = 3
YMAX_CHR_UNITS   = 25

ch_formatter = FuncFormatter(lambda v, pos: f"{v/CH_SCALE:.0f}")

AXES_LABEL_FONTSIZE = 18
TICK_LABEL_FONTSIZE = 14
TITLE_FONTSIZE      = 18
LEGEND_FONTSIZE     = 14

# ========= Smoothing Parameters =========
USE_RUNNING_AVG  = True
RA_WINDOW        = 12   # teacher requests w = 100 — adjust if desired

def apply_ch_axis(ax, ymax_units):
    ax.yaxis.set_major_formatter(ch_formatter)
    ax.set_yticks((np.arange(0, ymax_units + 1) * CH_SCALE).tolist())
    ax.set_ylim(0, ymax_units * CH_SCALE)
    ax.set_ylabel("CH(s) (×10⁻⁴)", fontsize=AXES_LABEL_FONTSIZE)
    ax.tick_params(axis='y', labelsize=TICK_LABEL_FONTSIZE)

def apply_x_axis(ax):
    ax.set_xlabel(f"Genomic distance ({X_UNIT})", fontsize=AXES_LABEL_FONTSIZE)
    ax.tick_params(axis='x', labelsize=TICK_LABEL_FONTSIZE)

def running_average(y, window):
    """
    Centered moving average ("teacher-required"):
      - np.convolve(..., mode='valid'): use only full windows
      - return length = len(y) - w + 1
    """
    y = np.asarray(y, dtype=float)
    L = len(y)
    if L == 0 or window <= 1:
        return y
    w = int(max(1, min(window, L)))
    kernel = np.ones(w, dtype=float) / w
    return np.convolve(y, kernel, mode='valid')

# ========= Helper Functions =========
def create_chromosome_labels(N, regions, labels):
    lab = [''] * N
    for (region, label) in zip(regions, labels):
        s, e = region
        for i in range(s, e):
            lab[i] = label
    return lab

def compute_row_normalized_P(matrix, chrom_labels, exclude_offset):
    """Mask near-diagonal interactions within same chromosome; keep inter-chrom; row-normalize."""
    N = matrix.shape[0]
    I, J = np.indices((N, N))
    chrom_arr   = np.array(chrom_labels)
    same_chrom  = (chrom_arr[:, None] == chrom_arr[None, :])
    intra_far   = np.abs(I - J) > exclude_offset
    mask_keep   = np.where(same_chrom, intra_far, True)

    A = matrix.astype(float).copy()
    A[~mask_keep] = 0.0
    row_sums = A.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return A / row_sums

def compute_Ps_both_triangles(P, N, exclude_offset):
    """
    Compute P(s) using both upper & lower triangle.
    Returns s, Ps(s), counts (= number of elements for each s).
    """
    s_vals = np.arange(exclude_offset + 1, N, dtype=int)
    Ps     = np.empty_like(s_vals, dtype=float); Ps.fill(np.nan)
    counts = np.zeros_like(s_vals, dtype=int)
    for idx, s in enumerate(s_vals):
        up   = np.diag(P,  k= s)
        low  = np.diag(P,  k=-s)
        vals = np.concatenate([up, low]) if (up.size or low.size) else np.array([])
        if vals.size > 0:
            Ps[idx]     = float(np.mean(vals))
            counts[idx] = vals.size
    valid = ~np.isnan(Ps)
    return s_vals[valid], Ps[valid], counts[valid]

def calc_CH_over_maps(P_list, N, exclude_offset):
    """
    Across multiple Hi-C maps:
      - compute Ps(s)
      - align to common length
      - CH(s) = standard deviation across maps
    """
    s_list, Ps_list = [], []
    for P in P_list:
        s, Ps, _ = compute_Ps_both_triangles(P, N, exclude_offset)
        if len(s) == 0:
            continue
        s_list.append(s)
        Ps_list.append(Ps)
    if not Ps_list:
        return np.array([]), np.array([])

    L = min(len(x) for x in s_list)
    s_common = s_list[0][:L]
    Ps_stack = np.vstack([arr[:L] for arr in Ps_list])  # (n_maps, L)
    CH = np.nanstd(Ps_stack, axis=0, ddof=1)
    return s_common, CH

# ========= Main Script =========
def main():
    chrom_labels_whole = create_chromosome_labels(N_BINS, REGIONS, LABELS)

    compare_whole = {}                      
    compare_by_chr = {lbl: {} for lbl in LABELS}

    for state, subdir in STATE_DIRS.items():
        folder = os.path.join(DATA_ROOT, subdir)
        files  = sorted(glob.glob(os.path.join(folder, GLOB_PATTERN)))
        if not files:
            print(f"[WARN] No files for state={state} in {folder}")
            continue
        print(f"[INFO] State={state}, files={len(files)}")

        P_list_whole = []
        P_lists_per_chr = {lbl: [] for lbl in LABELS}

        for fp in files:
            M = np.loadtxt(fp, delimiter=DELIM)

            # Whole genome
            P_whole = compute_row_normalized_P(M, chrom_labels_whole, EXCLUDE_OFFSET)
            P_list_whole.append(P_whole)

            # Per chromosome block
            for (start, end), lbl in zip(REGIONS, LABELS):
                block = M[start:end, start:end]
                block_labels = [lbl] * (end - start)
                P_block = compute_row_normalized_P(block, block_labels, EXCLUDE_OFFSET)
                P_lists_per_chr[lbl].append(P_block)

        # Whole-genome CH
        s_w, CH_w = calc_CH_over_maps(P_list_whole, N_BINS, EXCLUDE_OFFSET)
        if len(s_w):
            compare_whole[state] = (s_w, CH_w)

        # Per-chromosome CH
        for (start, end), lbl in zip(REGIONS, LABELS):
            L_chr = end - start
            if L_chr <= EXCLUDE_OFFSET + 1:
                continue
            s_c, CH_c = calc_CH_over_maps(P_lists_per_chr[lbl], L_chr, EXCLUDE_OFFSET)
            if len(s_c):
                compare_by_chr[lbl][state] = (s_c, CH_c)

    # ========= Plot WHOLE-genome CH(s) =========
    if compare_whole:
        min_len = min(len(v[0]) for v in compare_whole.values())
        fig, ax = plt.subplots(figsize=COMBINED_FIGSIZE)

        if USE_RUNNING_AVG:
            for state, (s, ch) in compare_whole.items():
                s_use  = s[:min_len]
                ch_use = ch[:min_len]
                ch_sm  = running_average(ch_use, RA_WINDOW)

                w_eff = int(max(1, min(RA_WINDOW, len(ch_use))))
                center_offset = (w_eff - 1) // 2
                s_center = s_use[center_offset : center_offset + len(ch_sm)]

                x_dist = bins_to_distance(s_center)
                label  = STATE_DISPLAY.get(state, state)
                ax.plot(x_dist, ch_sm, lw=1.8, marker='o', ms=3.0, label=label)
        else:
            for state, (s, ch) in compare_whole.items():
                s_use  = s[:min_len]
                x_dist = bins_to_distance(s_use)
                label  = STATE_DISPLAY.get(state, state)
                ax.plot(x_dist, ch[:min_len], lw=1.8, marker='o', ms=3.0, label=label)

        apply_x_axis(ax)
        apply_ch_axis(ax, YMAX_WHOLE_UNITS)

        ax.set_title("CH(s) Comparison Across States — Whole genome", fontsize=TITLE_FONTSIZE)
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)

        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "CH_comparison_WHOLE.png"),
                    dpi=300, bbox_inches="tight")
        plt.close(fig)

    # ========= Plot Per-Chromosome CH(s) =========
    for lbl in LABELS:
        states_for_lbl = compare_by_chr.get(lbl, {})
        if not states_for_lbl:
            continue
        min_len = min(len(v[0]) for v in states_for_lbl.values())
        fig, ax = plt.subplots(figsize=COMBINED_FIGSIZE)

        if USE_RUNNING_AVG:
            for state, (s, ch) in states_for_lbl.items():
                s_use  = s[:min_len]
                ch_use = ch[:min_len]
                ch_sm  = running_average(ch_use, RA_WINDOW)

                w_eff = int(max(1, min(RA_WINDOW, len(ch_use))))
                center_offset = (w_eff - 1) // 2
                s_center = s_use[center_offset : center_offset + len(ch_sm)]
                x_dist = bins_to_distance(s_center)

                label = STATE_DISPLAY.get(state, state)
                ax.plot(x_dist, ch_sm, lw=1.8, marker='o', ms=3.0, label=label)
        else:
            for state, (s, ch) in states_for_lbl.items():
                s_use  = s[:min_len]
                x_dist = bins_to_distance(s_use)
                label  = STATE_DISPLAY.get(state, state)
                ax.plot(x_dist, ch[:min_len], lw=1.8, marker='o', ms=3.0, label=label)

        apply_x_axis(ax)
        apply_ch_axis(ax, YMAX_CHR_UNITS)

        ax.set_title(f"CH(s) Comparison Across States — chr {lbl}",
                     fontsize=TITLE_FONTSIZE)
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE)

        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"CH_comparison_chr_{lbl}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(f"[DONE] Comparison plots saved in: {OUT_DIR}")
    print(f"      Running average: {'ON' if USE_RUNNING_AVG else 'OFF'}, window={RA_WINDOW}")
    print(f"      X-axis unit: {X_UNIT}, BIN_SIZE_BP={BIN_SIZE_BP}")

if __name__ == "__main__":
    main()
