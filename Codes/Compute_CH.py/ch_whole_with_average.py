# -*- coding: utf-8 -*-
"""
CH(s) from Hi-C probabilities using row-normalized P.
Outputs ONLY:
  - CH_comparison_WHOLE.png (whole-genome, across states)
  - CH_comparison_chr_<LABEL>.png for LABEL in [2L,2R,3L,3R,4,X] (across states)

Rules:
  - Teacher-approved row normalization: mask same-chrom |i-j| <= EXCLUDE_OFFSET; keep inter-chrom
  - Use BOTH upper & lower triangles when computing P(s)
  - Y axis displayed as 0–8 with label "CH(s) (×10⁻⁴)"
  - Running average smoothing (window=RA_WINDOW, 10–50), edge-padded, same length

Author: You
Date: 2025-10-30
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ========= CONFIGURATION =========
DATA_ROOT = r"D:\Chro"                 # ← Change to your data root
OUT_DIR   = r"./out_ch_running_average_4l"   # Output directory
os.makedirs(OUT_DIR, exist_ok=True)

STATE_DIRS = {
    "WT_orig":        "WT3h_orig",
    "CD3h":           "CD3h",
    "LM3h":           "LM3h",
    "WT_after_CD3h":  "WT3h_after_CD3h",
    "WT_after_LM3h":  "WT3h_after_LM3h",
}
GLOB_PATTERN = "*.map"
DELIM = ','

EXCLUDE_OFFSET = 4                     # Near-diagonal exclusion width
N_BINS  = 1169                         # Total genome bins
REGIONS = [(0, 225), (225, 438), (438, 659), (659, 966), (966, 985), (985, N_BINS)]
LABELS  = ["2L", "2R", "3L", "3R", "4", "X"]

# Plotting and y-axis settings
plt.rcParams["figure.facecolor"] = "white"
COMBINED_FIGSIZE = (7.5, 5)
CH_SCALE         = 1e-4               # Display 0–8 × 10⁻⁴ as 0–8
YMAX_UNITS       = 6
ch_formatter     = FuncFormatter(lambda v, pos: f"{v/CH_SCALE:.0f}")

# Smoothing parameters (recommended 10–50)
USE_RUNNING_AVG  = True
RA_WINDOW        = 50                  # ← Adjust between 10–50; disable smoothing by setting USE_RUNNING_AVG = False above

def apply_ch_axis(ax):
    ax.yaxis.set_major_formatter(ch_formatter)
    ax.set_yticks((np.arange(0, YMAX_UNITS + 1) * CH_SCALE).tolist())
    ax.set_ylim(0, YMAX_UNITS * CH_SCALE)
    ax.set_ylabel("CH(s) (×10⁻⁴)")

def running_average(y, window):
    """Running average with edge padding, output has same length."""
    y = np.asarray(y, dtype=float)
    L = len(y)
    if L == 0 or window <= 1:
        return y
    w = int(max(1, min(window, L)))      # clamp to [1, L]
    kernel = np.ones(w, dtype=float) / w
    pad_left = w // 2
    pad_right = w - 1 - pad_left
    ypad = np.pad(y, (pad_left, pad_right), mode='edge')
    return np.convolve(ypad, kernel, mode='valid')

# ========= Utility functions =========
def create_chromosome_labels(N, regions, labels):
    """Assign chromosome labels for each genome bin."""
    lab = [''] * N
    for (region, label) in zip(regions, labels):
        s, e = region
        for i in range(s, e):
            lab[i] = label
    return lab

def compute_row_normalized_P(matrix, chrom_labels, exclude_offset):
    """
    Row-normalize with masking:
      - same chromosome & |i-j| <= exclude_offset: masked
      - inter-chromosomal interactions kept
    """
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
    Compute P(s) using both upper and lower triangles:
      s = exclude_offset+1 .. N-1
    Returns: s_values, Ps(s), counts = 2*(N-s)
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
    For multiple maps:
      - compute Ps(s) for each map
      - align to common min length
      - CH(s) = std over maps of Ps(s)
    Returns s_common, CH(s)
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

# ========= Main pipeline =========
def main():
    chrom_labels_whole = create_chromosome_labels(N_BINS, REGIONS, LABELS)

    # Across-state comparison containers
    compare_whole = {}                           # state -> (s, CH)
    compare_by_chr = {lbl: {} for lbl in LABELS} # chr -> {state -> (s, CH)}

    for state, subdir in STATE_DIRS.items():
        folder = os.path.join(DATA_ROOT, subdir)
        files  = sorted(glob.glob(os.path.join(folder, GLOB_PATTERN)))
        if not files:
            print(f"[WARN] No files for state={state} in {folder}")
            continue
        print(f"[INFO] State={state}, files={len(files)}")

        # Collect P for whole genome & per chromosome
        P_list_whole = []
        P_lists_per_chr = {lbl: [] for lbl in LABELS}

        for fp in files:
            M = np.loadtxt(fp, delimiter=DELIM)

            # Whole genome
            P_whole = compute_row_normalized_P(M, chrom_labels_whole, EXCLUDE_OFFSET)
            P_list_whole.append(P_whole)

            # Chromosome blocks
            for (start, end), lbl in zip(REGIONS, LABELS):
                block = M[start:end, start:end]
                block_labels = [lbl] * (end - start)
                P_block = compute_row_normalized_P(block, block_labels, EXCLUDE_OFFSET)
                P_lists_per_chr[lbl].append(P_block)

        # Whole-genome CH(s)
        s_w, CH_w = calc_CH_over_maps(P_list_whole, N_BINS, EXCLUDE_OFFSET)
        if len(s_w):
            compare_whole[state] = (s_w, CH_w)

        # Per-chromosome CH(s)
        for (start, end), lbl in zip(REGIONS, LABELS):
            L = end - start
            if L <= EXCLUDE_OFFSET + 1:
                continue
            s_c, CH_c = calc_CH_over_maps(P_lists_per_chr[lbl], L, EXCLUDE_OFFSET)
            if len(s_c):
                compare_by_chr[lbl][state] = (s_c, CH_c)

    # —— Plot WHOLE (across states) —— #
    if compare_whole:
        min_len = min(len(v[0]) for v in compare_whole.values())
        fig, ax = plt.subplots(figsize=COMBINED_FIGSIZE)
        for state, (s, ch) in compare_whole.items():
            ch_use = ch[:min_len]
            if USE_RUNNING_AVG:
                ch_use = running_average(ch_use, RA_WINDOW)
            ax.plot(s[:min_len], ch_use, lw=1.8, marker='o', ms=3.0, label=state)
        ax.set_xlabel("Genomic distance s (bins)")
        apply_ch_axis(ax)
        title = "CH(s) Comparison Across States — WHOLE (both triangles)"
        if USE_RUNNING_AVG:
            title += f" | running avg (w={RA_WINDOW})"
        ax.set_title(title)
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "CH_comparison_WHOLE.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    # —— Plot per-chromosome comparisons —— #
    for lbl in LABELS:
        states_for_lbl = compare_by_chr.get(lbl, {})
        if not states_for_lbl:
            continue
        min_len = min(len(v[0]) for v in states_for_lbl.values())
        fig, ax = plt.subplots(figsize=COMBINED_FIGSIZE)
        for state, (s, ch) in states_for_lbl.items():
            ch_use = ch[:min_len]
            if USE_RUNNING_AVG:
                ch_use = running_average(ch_use, RA_WINDOW)
            ax.plot(s[:min_len], ch_use, lw=1.8, marker='o', ms=3.0, label=state)
        ax.set_xlabel("Genomic distance s (bins)")
        apply_ch_axis(ax)
        title = f"CH(s) Comparison Across States — {lbl} (both triangles)"
        if USE_RUNNING_AVG:
            title += f" | running avg (w={RA_WINDOW})"
        ax.set_title(title)
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"CH_comparison_chr_{lbl}.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(f"[DONE] Comparison plots saved in: {OUT_DIR}")
    print(f"      Running average: {'ON' if USE_RUNNING_AVG else 'OFF'}, window={RA_WINDOW}")

if __name__ == "__main__":
    main()
