# -*- coding: utf-8 -*-
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# =========================
# 可调参数（统一样式）
# =========================
OUT_DIR = r"D:\Chro"
os.makedirs(OUT_DIR, exist_ok=True)

# === 通用样式设置 ===
BASE_FS   = float(plt.rcParams.get('font.size', 10.0))
LABEL_FS  = BASE_FS * 3
TICK_FS   = BASE_FS * 2
LEGEND_FS = BASE_FS * 1.2

MAJOR_LEN = float(plt.rcParams.get('xtick.major.size', 3.5)) * 2
MAJOR_W   = float(plt.rcParams.get('xtick.major.width', 0.8)) * 2
MINOR_LEN = float(plt.rcParams.get('xtick.minor.size', 2.0)) * 2
MINOR_W   = float(plt.rcParams.get('xtick.minor.width', 0.6)) * 2

transition_color = 'C3'
connect_style = dict(linestyle='--', linewidth=1.5, alpha=0.9, color=transition_color)

k_formatter = FuncFormatter(lambda v, pos: f"{v/1000:.1f}")

def apply_common_axis_style(ax, ylim=(9000, 10600), y_label='Entropy (×10^3 bits)'):
    ax.set_ylim(*ylim)
    ax.yaxis.set_major_formatter(k_formatter)
    ax.set_xlabel('Time (minutes)', fontsize=LABEL_FS)
    ax.set_ylabel(y_label, fontsize=LABEL_FS)
    ax.tick_params(axis='both', which='major',
                   labelsize=TICK_FS, length=MAJOR_LEN, width=MAJOR_W)
    ax.tick_params(axis='both', which='minor',
                   labelsize=TICK_FS, length=MINOR_LEN, width=MINOR_W)

# =========================
# 数据：WT → Lamins depleted → WT (Rejuvenated)
# =========================

# WT (original)
wt_min = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
wt_total = [
    9084.4452,9062.0980,9146.3733,9143.9850,9067.5422,9088.7596,9129.6872,9076.6409,
    9121.8193,9121.5474,9074.9941,9147.1354,9089.7314,9113.3103,9109.8950,9077.1823,
    9120.0304,9123.5952
]

# WT → Lamins depleted (transition)
trans_wt_lm_min = list(range(1,13))
trans_wt_lm_total = [
    9835.5433,9900.7579,9927.9688,9932.2644,9912.0372,9910.8702,
    9907.7002,9901.1786,9910.0985,9912.0000,9947.9202,9914.0747
]

# Lamins depleted
lm_min = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
lm_total = [
    9924.8077,9905.5001,9960.7415,9944.3725,10003.5202,9962.5103,9944.4652,9947.3499,
    9901.7054,9964.3221,9904.1504,9948.0381,9940.0597,9905.5784,9936.5336,9958.6023,
    9948.8370,9941.3637
]

# Lamins depleted → Rejuvenated WT (transition)
trans_lm_rejuv_min = list(range(1,13))
trans_lm_rejuv_total = [
    9426.5923,9216.4047,9168.5095,9082.5544,9085.7436,9103.4357,
    9106.4794,9112.1292,9094.2405,9071.8116,9137.5752,9128.9585
]

# Rejuvenated WT
rejuv_min = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
rejuv_total = [
    9072.3546,9109.3698,9125.1919,9094.7952,9150.8516,9088.4210,9051.8155,9119.1590,
    9111.8105,9052.9421,9112.6178,9132.2148,9087.0136,9059.2901,9156.0340,9119.9448,
    9121.7281,9124.0032
]

# =========================
# 串联时间轴
# =========================
wt_x       = wt_min
trans1_x   = [m + 180 for m in trans_wt_lm_min]
lm_x       = [m + 192 for m in lm_min]
trans2_x   = [m + 372 for m in trans_lm_rejuv_min]
rejuv_x    = [m + 384 for m in rejuv_min]

# =========================
# 绘图
# =========================
plt.figure(figsize=(12,6))

# 主阶段曲线
plt.plot(wt_x,    wt_total,     label='WT (Original)',           marker='o', linestyle='-')
plt.plot(lm_x,    lm_total,     label='Lamins depleted',         marker='o', linestyle='-')
plt.plot(rejuv_x, rejuv_total,  label='WT (Rejuvenated)',        marker='o', linestyle='-')

# 两段过渡（均进图例）
plt.plot(trans1_x, trans_wt_lm_total, linestyle='--', marker=None, color=transition_color,
         label='Transition (WT → Lamins depleted)')
plt.plot(trans2_x, trans_lm_rejuv_total, linestyle='--', marker=None, color=transition_color,
         label='Transition (Lamins depleted → Rejuvenated WT)')

# 竖直连接（也为虚线）
plt.plot([wt_x[-1],     trans1_x[0]], [wt_total[-1],     trans_wt_lm_total[0]],  **connect_style)
plt.plot([trans1_x[-1], lm_x[0]],     [trans_wt_lm_total[-1], lm_total[0]],      **connect_style)
plt.plot([lm_x[-1],     trans2_x[0]], [lm_total[-1],     trans_lm_rejuv_total[0]],  **connect_style)
plt.plot([trans2_x[-1], rejuv_x[0]],  [trans_lm_rejuv_total[-1], rejuv_total[0]],   **connect_style)

# 过渡散点
plt.scatter(trans1_x, trans_wt_lm_total, s=12, marker='o', color=transition_color, zorder=5)
plt.scatter(trans2_x, trans_lm_rejuv_total, s=12, marker='o', color=transition_color, zorder=5)

# 标题
plt.title('WT → Lamins depleted → WT (Rejuvenated): Entropy')

# 统一样式
ax = plt.gca()
apply_common_axis_style(ax)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=LEGEND_FS)

# 保存
out_path = os.path.join(OUT_DIR, "total_entropy_WT_LaminsDepleted_RejuvenatedWT.png")
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close()

print("Saved:", out_path)
