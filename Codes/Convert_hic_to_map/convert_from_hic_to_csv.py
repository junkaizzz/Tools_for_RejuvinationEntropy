"""
Author: Junkai Zhang (Kain)
Created:  '2025-09-15'
Last Modified:2025-09-15
Email: junkaiz@vt.edu
"""
import numpy as np
import pandas as pd

# Step 1: Read raw data
df = pd.read_csv("mouse_cluster5_chr11_250kb.txt", sep="\t", header=None) # Use our own data
bin_size = 250000

# Step 2: Calculate bin indices
df["i"] = df[0] // bin_size
df["j"] = df[1] // bin_size

# Step 3: Create an empty matrix
size = max(df["i"].max(), df["j"].max()) + 1
matrix = np.zeros((size, size))

# Step 4: Fill the matrix (symmetric)
for _, row in df.iterrows():
    i = int(row["i"])
    j = int(row["j"])
    value = row[2]
    matrix[i][j] = value
    matrix[j][i] = value

# Step 5: Save as CSV
np.savetxt("mouse_cluster5_chr11_250kb.csv", matrix, delimiter=",")

# (Optional) Save as .npy format (faster)
# np.save("output_chr1_matrix.npy", matrix)

print("complete")
print("Matrix size:", matrix.shape)
