import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Path to results
folder = '/Users/axdenis03/Documents/University /BioInformatics/Project_L11/results'

# Get all files
all_files = os.listdir(folder)

# Extract only the original fasta names from the _vs_ files
sequence_files = set()
for f in all_files:
    if "_vs_" in f and f.endswith('.txt'):
        seq1, seq2 = f.split("_vs_")
        seq2 = seq2.replace(".txt","")
        sequence_files.add(seq1)
        sequence_files.add(seq2)

sequence_files = sorted(sequence_files)
print("Sequences found:", sequence_files)
n = len(sequence_files)

# Map file names to indices
file_indices = {name: idx for idx, name in enumerate(sequence_files)}

# Initialize similarity matrix
similarity_matrix = np.zeros((n, n))
for i in range(n):
    similarity_matrix[i, i] = 100

# Fill matrix from results files
for f in all_files:
    if f.endswith('.txt') and "_vs_" in f:
        seq1, seq2 = f.split("_vs_")
        seq2 = seq2.replace(".txt","")
        idx1, idx2 = file_indices[seq1], file_indices[seq2]
        with open(os.path.join(folder,f), 'r') as file:
            content = file.read()
            match = re.search(r"Similarity:\s*([\d\.]+)%", content)
            if match:
                sim = float(match.group(1))
                similarity_matrix[idx1, idx2] = sim
                similarity_matrix[idx2, idx1] = sim
            else:
                print(f"Warning: no similarity found in {f}")

# Convert to DataFrame
df = pd.DataFrame(similarity_matrix, index=sequence_files, columns=sequence_files)
print(df)

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(df, annot=True, fmt=".2f", cmap="viridis", cbar_kws={'label': 'Similarity (%)'})
plt.title("Pairwise Genome Similarity Heatmap")
plt.xticks(rotation=45)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()