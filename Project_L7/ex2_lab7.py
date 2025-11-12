# 4. Download 10 influenza genomes. For each, plot the frequency of consecutive repeats of length 6-10.

import os
from collections import Counter
import matplotlib.pyplot as plt

def read_fasta(file_path):
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()

def find_consecutive_repeats(seq, min_len=6, max_len=10):
    repeats = {}
    for k in range(min_len, max_len + 1):
        counts = Counter()
        i = 0
        while i <= len(seq) - k:
            kmer = seq[i:i+k]
            repeat_count = 0
            while seq[i + repeat_count*k : i + (repeat_count+1)*k] == kmer:
                repeat_count += 1
            if repeat_count >= 2:
                counts[kmer] += 1  # just count it once per consecutive occurrence
                i += repeat_count * k
            else:
                i += 1
        repeats[k] = counts
    return repeats

folder = '/Users/axdenis03/Documents/University /BioInformatics/Project_L7'

fig, axes = plt.subplots(2, 5, figsize=(25, 10))
axes = axes.flatten()

for n in range(1, 11):
    file_path = os.path.join(folder, f"sequence-{n}.fasta")
    seq = read_fasta(file_path)
    repeats = find_consecutive_repeats(seq, 6, 10)

    x_labels = []
    y_counts = []
    for k in sorted(repeats.keys()):
        counter = repeats[k]
        for motif, count in counter.items():
            x_labels.append(f"{motif} ({k}-bp)")
            y_counts.append(count)

    ax = axes[n-1]
    if x_labels:
        ax.bar(range(len(y_counts)), y_counts, color="#1f77b4")
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=90, fontsize=7)
        ax.set_ylabel("Occurrences")
        ax.set_title(f"Genome {n}")
    else:
        ax.set_title(f"Genome {n} (no repeats)")

plt.tight_layout()
plt.show()
