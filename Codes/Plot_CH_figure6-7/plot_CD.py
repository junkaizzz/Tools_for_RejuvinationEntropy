# -*- coding: utf-8 -*-
"""
Author: Junkai Zhang (Kain)
Created:  '2025-09-15'
Last Modified: 2025-09-15
Email: junkaiz@vt.edu
"""

import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# =========================
# Adjustable parameters (Unified style)
# =========================
OUT_DIR = r"D:\Chro"
os.makedirs(OUT_DIR, exist_ok=True)

# Unified font sizes (using numeric base value to avoid rcParams returning strings)
BASE_FS   = float(plt.rcParams.get('font.size', 10.0))   # Base font size
LABEL_FS  = BASE_FS * 3                                  # Axis label ×3
TICK_FS   = BASE_FS * 2                                  # Tick label ×2
LEGEND_FS = BASE_FS * 1.2                                # Legend font size

# Major/minor tick length & width doubled
MAJOR_LEN = float(plt.rcParams.get('xtick.major.size', 3.5)) * 2
MAJOR_W   = float(plt.rcParams.get('xtick.major.width', 0.8)) * 2
MINOR_LEN = float(plt.rcParams.get('xtick.minor.size', 2.0)) * 2
MINOR_W   = float(plt.rcParams.get('xtick.minor.width', 0.6)) * 2

# Transition segments and vertical connectors: unified dashed line style
transition_color = 'C3'
connect_style = dict(linestyle='--', linewidth=1.5, alpha=0.9, color=transition_color)

# y-axis formatter for 10^3 bits (e.g., 9.0, 9.1, ...)
k_formatter = FuncFormatter(lambda v, pos: f"{v/1000:.1f}")

def apply_common_axis_style(ax, ylim=(9000, 10600), y_label='Entropy (×10^3 bits)'):
    """Apply shared axis formatting for both systems."""
    ax.set_ylim(*ylim)                                     # Fixed y-axis range
    ax.yaxis.set_major_formatter(k_formatter)              # Format y-axis as 10^3 bits
    ax.set_xlabel('Time (minutes)', fontsize=LABEL_FS)
    ax.set_ylabel(y_label, fontsize=LABEL_FS)

    # Tick labels and tick mark sizes ×2
    ax.tick_params(axis='both', which='major',
                   labelsize=TICK_FS, length=MAJOR_LEN, width=MAJOR_W)
    ax.tick_params(axis='both', which='minor',
                   labelsize=TICK_FS, length=MINOR_LEN, width=MINOR_W)

# =========================
# Data (Total entropy)
# =========================
# WT (original)
wt_min = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
wt_total = [9084.4452,9062.0980,9146.3733,9143.9850,9067.5422,9088.7596,9129.6872,9076.6409,
            9121.8193,9121.5474,9074.9941,9147.1354,9089.7314,9113.3103,9109.8950,9077.1823,
            9120.0304,9123.5952]

# WT → Dead transition (1–12 min)
trans_wt_dead_min = list(range(1,13))
trans_wt_dead_total = [10319.4303,10405.6616,10449.7950,10462.4452,10451.3524,10450.8085,
                       10470.4791,10462.9176,10466.8809,10467.5281,10472.8894,10456.4492]

# Dead (CD)
cd_min = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
cd_total = [10467.7340,10472.6168,10464.5179,10449.1549,10459.7882,10452.9637,10478.3798,10422.8005,
            10459.9942,10487.4076,10479.6606,10457.8396,10423.7628,10432.8876,10454.7405,10469.9840,
            10503.0664,10470.2304]

# Dead → Rejuvenated WT transition (1–12 min)
trans_rejuv_min = list(range(1,13))
trans_rejuv_total = [9467.1823,9186.6456,9179.9947,9138.4908,9144.8765,9144.8359,9147.8627,9075.4063,
                     9104.7279,9117.0151,9106.7121,9179.5841]

# WT (rejuvenated)
wt_after_min = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
wt_after_total = [9138.9970,9082.4718,9113.2729,9151.0673,9100.8204,9097.8799,9043.4454,9114.7609,
                  9096.3981,9131.2399,9144.6757,9131.4584,9152.9074,9123.7727,9105.1475,9091.2796,
                  9126.2022,9100.3961]

# =========================
# Serial time axis (same logic as before)
# =========================
wt_x       = wt_min
trans1_x   = [m + 180 for m in trans_wt_dead_min]   # 181..192
cd_x       = [m + 192 for m in cd_min]              # 202..372
trans2_x   = [m + 372 for m in trans_rejuv_min]     # 373..384
wt_after_x = [m + 384 for m in wt_after_min]        # 394..564

# =========================
# Plotting
# =========================
plt.figure(figsize=(12,6))

# Main curves
l_wt,      = plt.plot(wt_x,      wt_total,            label='WT (Original)',       marker='o', linestyle='-')
l_cd,      = plt.plot(cd_x,      cd_total,            label='CD (Dead)',           marker='o', linestyle='-')
l_after,   = plt.plot(wt_after_x,wt_after_total,      label='WT (Rejuvenated)',    marker='o', linestyle='-')

# Transition segments (included in legend)
l_trans1,  = plt.plot(trans1_x,  trans_wt_dead_total, label='Transition (WT → Dead)',
                      linestyle='--', marker=None, color=transition_color)
l_trans2,  = plt.plot(trans2_x,  trans_rejuv_total,   label='Transition (Dead → Rejuvenated WT)',
                      linestyle='--', marker=None, color=transition_color)

# Vertical dashed connectors
plt.plot([wt_x[-1],     trans1_x[0]], [wt_total[-1],     trans_wt_dead_total[0]], **connect_style)
plt.plot([trans1_x[-1], cd_x[0]],     [trans_wt_dead_total[-1], cd_total[0]],     **connect_style)
plt.plot([cd_x[-1],     trans2_x[0]], [cd_total[-1],     trans_rejuv_total[0]],   **connect_style)
plt.plot([trans2_x[-1], wt_after_x[0]],[trans_rejuv_total[-1], wt_after_total[0]],**connect_style)

# Transition markers
plt.scatter(trans1_x, trans_wt_dead_total, s=10, marker='o', color=l_trans1.get_color(), zorder=5)
plt.scatter(trans2_x, trans_rejuv_total,   s=10, marker='o', color=l_trans2.get_color(), zorder=5)

# Title
plt.title('WT → Dead → WT (Rejuvenated): Entropy')

# Axis & legend formatting
ax = plt.gca()
apply_common_axis_style(ax)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=LEGEND_FS)

# Save figure
out_path = os.path.join(OUT_DIR, "total_entropy_all_points.png")
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close()

print("Saved:", out_path)
