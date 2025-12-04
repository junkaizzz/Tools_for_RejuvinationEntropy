# Toolkit

Scripts for converting between `.map` and `.hic` formats and combining two Hi-C maps, built around **Juicer Tools**, plus new analytical scripts for **row normalization**, **entropy**, and **CH(s)** computation.

---

## Requirements

- **Java JDK 1.7 or 1.8**
  Download: <http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html>
  (Ubuntu/LinuxMint alternative: <http://tecadmin.net/install-oracle-java-8-jdk-8-ubuntu-via-ppa/>)
  Minimum system requirements for Java: <http://java.com/en/download/help/sysreq.xml>

- **Juicer Tools (latest .jar)**
  <https://github.com/aidenlab/juicer/wiki/Download>

- **Python 3.9+**
  Dependencies: `numpy`, `pathlib`, `matplotlib`
  ```bash
  pip install numpy matplotlib pathlib
  ```

---

## Usage

### 1) `.map → .hic`

1. Edit `Convert_map_to_hic/convert_map_to_hic.py` to specify the `.map` files.
2. Run the script to produce the `.hic` file.

### 2) `.hic → .csv (.map)`

1. Use Juicer Tools to dump the observed matrix, e.g.:
   ```bash
   java -jar juicer_tools.jar dump observed NONE input.hic chr1 chr1 BP 250000 output.txt
   ```
2. Run `convert_from_hic_to_csv.py` to convert to `.csv` (or `.map`).

### 3) Combine Two `.hic` Maps

Use `Combine_2hic.py` to merge two `.hic` files into one composite view.

---

## 4) Row-normalization & Entropy

**Script:** `row_norm_entropy.py`

Performs **row normalization** on Hi-C matrices (masking near-diagonal intra-chromosomal regions, preserving inter-chromosomal interactions) and computes **entropy** for each row and globally.

### Inputs
- `.map` matrix (comma-delimited)
- Chromosome region definitions (`regions`) and labels (`labels`)

### Parameters
- `exclude_offset`: width of diagonal mask (default = 4)
- `file_path`: path to the input `.map`

### Outputs
- Console: total entropy and entropy per TAD (bits)
- Figure: `hic_matrix_heatmap.png` — log-scale P-matrix heatmap

### Example
```bash
python row_norm_entropy.py
```

---

## 5) CH(s) Comparison Across States

**Script:** `ch_whole_with_average.py`

Computes **CH(s)** (Contact Heterogeneity) across multiple biological states using **row-normalized P matrices**.

### Directory Structure Example
```
DATA_ROOT/
  WT3h_orig/
  CD3h/
  LM3h/
  WT3h_after_CD3h/
  WT3h_after_LM3h/
```

### Configuration Highlights
- `EXCLUDE_OFFSET = 4`
- `USE_RUNNING_AVG = True`
- `RA_WINDOW = 50`
- `REGIONS` and `LABELS` match your genome segmentation

### Outputs
- `CH_comparison_WHOLE.png` — whole-genome comparison
- `CH_comparison_chr_<LABEL>.png` — per-chromosome comparison

### Example
```bash
python ch_whole_with_average.py
```

---

## Reference

- Juicer Tools Documentation: <https://github.com/aidenlab/juicer/wiki/Pre>

