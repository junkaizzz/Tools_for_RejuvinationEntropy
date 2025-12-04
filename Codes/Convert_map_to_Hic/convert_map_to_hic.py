"""
Author: Junkai Zhang (Kain)
Created:  '2025-09-15'
Last Modified:2025-09-15
Email: junkaiz@vt.edu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert Igor “.map” matrices ➜ Juicer short‑with‑score (TAB) ➜ “.hic”.
Requirements: numpy, pandas, gzip (Python) + Java + juicer_tools_1.22.01.jar
"""

# ======== CONFIGURATION (edit only this block) ========
BIN_SIZE = 10_000          # Resolution of Igor matrix in base pairs
SCALE    = 1e6             # If values are frequencies (0‑1), multiply by SCALE to make integers;
                           # if values are already integer counts, set SCALE = 1

MAP_FILES = [              # One or more .map files to process
    r"Week7\Hic-map_new\WT_1min_orig_HiC_0.2_av18.map"
]

# Bin ranges (left‑closed, right‑open) and their corresponding chromosome names
chrom_regions_list = [(0,225), (225,438), (438,659), (659,966), (966,985), (985,1169)]
chrom_labels_names = ["2L", "2R", "3L", "3R", "4", "X"]
# =====================================================

import numpy as np, pandas as pd, gzip, subprocess, sys
from pathlib import Path

JAR = Path("Tools/juicer_tools_1.22.01.jar")   # Path to juicer_tools JAR
RAM = "4g"                                     # Java heap size

# ---------- build lookup table: bin → (chr, start) ----------
def build_lookup():
    rows = []
    for (b0, b1), chrom in zip(chrom_regions_list, chrom_labels_names):
        for b in range(b0, b1):
            rows.append((b, chrom, (b - b0) * BIN_SIZE))
    return pd.DataFrame(rows, columns=["bin", "chr", "start"])

LOOKUP = build_lookup()

# ---------- process a single .map file ----------
def process_one(fname):
    m = Path(fname)
    if not m.exists():
        print(f"[WARN] file not found: {m}")
        return
    print(f"\n=== processing {m.name} ===")

    # load matrix, keep upper triangle (including diagonal)
    mat = np.loadtxt(m, delimiter=",")
    r, c = np.triu_indices_from(mat, k=0)
    v    = mat[r, c]
    keep = v > 0
    r, c, v = r[keep], c[keep], v[keep]

    # scale and cast to 32‑bit int
    if SCALE != 1:
        v = np.rint(v * SCALE)
    v = v.astype(np.int32)

    # merge coordinates
    df = pd.DataFrame({"bin1": r, "bin2": c, "count": v})
    df = df.merge(LOOKUP, left_on="bin1", right_on="bin") \
           .rename(columns={"chr": "chr1", "start": "pos1"}).drop(columns="bin")
    df = df.merge(LOOKUP, left_on="bin2", right_on="bin") \
           .rename(columns={"chr": "chr2", "start": "pos2"}).drop(columns="bin")

    # Juicer requires rows sorted by chr1 → chr2 → pos1 → pos2
    df = df.sort_values(["chr1", "chr2", "pos1", "pos2"])

    # write 9‑column short‑with‑score file
    out_txt = m.with_suffix(".short_score.txt.gz")
    with gzip.open(out_txt, "wt") as gz:
        for row in df.itertuples(index=False):
            # str1 chr1 pos1 frag1  str2 chr2 pos2 frag2  score
            gz.write(
                f"0\t{row.chr1}\t{row.pos1}\t0\t"
                f"0\t{row.chr2}\t{row.pos2}\t1\t"
                f"{row.count}\n"
            )
    print(f"  · wrote {len(df):,} rows ➜ {out_txt}")

    # write chrom.sizes
    sizes_file = m.with_suffix(".chrom.sizes")
    with open(sizes_file, "w") as f:
        for (b0, b1), chrom in zip(chrom_regions_list, chrom_labels_names):
            f.write(f"{chrom}\t{(b1 - b0) * BIN_SIZE}\n")
    print(f"  · chrom.sizes ➜ {sizes_file}")

    # run Juicer pre
    hic_out = m.with_suffix(".hic")
    cmd = ["java", f"-Xmx{RAM}", "-jar", str(JAR),
           "pre", "-r", str(BIN_SIZE),
           str(out_txt), str(hic_out), str(sizes_file)]
    print("  · juicer:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"=== DONE → {hic_out}")

# ---------- entry point ----------
def main():
    for f in MAP_FILES:
        process_one(f)

if __name__ == "__main__":
    main()
