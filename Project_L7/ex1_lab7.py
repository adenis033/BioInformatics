# 1. Load a DNA sequence (1000–3000 bases) from FASTA.
# 2. Detect consecutive repeated motifs of length 6–10.
# 3. Plot one bar chart per motif length showing motifs vs. number of consecutive repeats.

import random
from collections import Counter
import matplotlib.pyplot as plt

def load_fasta(path):
    with open(path, "r") as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()

def extract_random_region(seq, min_len=1000, max_len=3000):
    if len(seq) <= max_len:
        return seq
    start = random.randint(0, len(seq) - max_len)
    length = random.randint(min_len, max_len)
    return seq[start:start+length]

def find_consecutive_repeats(seq, min_k=6, max_k=10):
    results = {}
    for k in range(min_k, max_k+1):
        counter = Counter()
        i = 0
        while i <= len(seq) - k:
            motif = seq[i:i+k]
            repeats = 1
            j = i + k
            while j <= len(seq) - k and seq[j:j+k] == motif:
                repeats += 1
                j += k
            if repeats > 1:
                counter[motif] = max(counter.get(motif, 0), repeats)
            i += 1
        results[k] = counter
    return results

def plot_repeats_all(repeats_dict):
    combined = {}
    for k in repeats_dict:
        combined.update(repeats_dict[k])
    
    if not combined:
        print("No consecutive repeats found.")
        return
    
    motifs, counts = zip(*sorted(combined.items(), key=lambda x: x[1], reverse=True))
    
    plt.figure(figsize=(14,6))
    bars = plt.bar(range(len(motifs)), counts, color="#4B9CD3", edgecolor="black")
    
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, count + 0.1, str(count),
                 ha='center', va='bottom', fontsize=8)
    
    plt.xticks(range(len(motifs)), motifs, rotation=90, fontsize=8)
    plt.ylabel("Consecutive repeats")
    plt.title("Consecutive DNA repeats (6–10 bases)")
    plt.yticks(range(1, max(counts)+2))  # integer y-axis
    plt.tight_layout()
    plt.show()

# main execution
fasta_file = '/Users/axdenis03/Documents/University /BioInformatics/Project_L7/sequence.fasta'
full_seq = load_fasta(fasta_file)
dna_seq = extract_random_region(full_seq)

repeats = find_consecutive_repeats(dna_seq, 6, 10)
plot_repeats_all(repeats)
