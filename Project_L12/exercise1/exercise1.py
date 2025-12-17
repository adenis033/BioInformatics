import numpy as np
import matplotlib.pyplot as plt

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT"
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"
bases = ['A', 'C', 'G', 'T']
motif_len = len(motifs[0])
num_seqs = len(motifs)

count_matrix = np.zeros((4, motif_len))
freq_matrix = np.zeros((4, motif_len))
log_lik_matrix = np.zeros((4, motif_len))

for seq in motifs:
    for i, nucleotide in enumerate(seq):
        idx = bases.index(nucleotide)
        count_matrix[idx, i] += 1

for r in range(4):
    for c in range(motif_len):
        freq = count_matrix[r, c] / num_seqs
        freq_matrix[r, c] = freq
        
        p_n = (count_matrix[r, c] + 0.01) / (num_seqs + 0.04) 
        log_lik_matrix[r, c] = np.log(p_n / 0.25)


scan_scores = []
positions = []

for i in range(len(S) - motif_len + 1):
    window = S[i : i + motif_len]
    score = 0
    for pos, nuc in enumerate(window):
        if nuc in bases:
            row = bases.index(nuc)
            score += log_lik_matrix[row, pos]
    
    scan_scores.append(score)
    positions.append(i)

print("--- Count Matrix ---")
print("   " + "  ".join([str(i+1) for i in range(motif_len)]))
for i, base in enumerate(bases):
    print(f"{base}: {count_matrix[i]}")

print("\n--- Relative Frequencies Matrix ---")
for i, base in enumerate(bases):
    print(f"{base}: {np.round(freq_matrix[i], 2)}")

print("\n--- Log-Likelihood Matrix ---")
for i, base in enumerate(bases):
    print(f"{base}: {np.round(log_lik_matrix[i], 2)}")

print("\n--- Max Score in Sequence S ---")
max_score = max(scan_scores)
max_pos = np.argmax(scan_scores)
print(f"Max Signal: {max_score:.4f} at Index {max_pos}")
print(f"Subsequence: {S[max_pos : max_pos + motif_len]}")

plt.figure(figsize=(10, 6))
plt.plot(positions, scan_scores, marker='o', color='b', linewidth=1.5, label='Log-Likelihood Score')
plt.axhline(0, color='red', linestyle='--', linewidth=1, label='Neutral Threshold')
plt.xlabel('Position in Sequence S')
plt.ylabel('Score')
plt.title('Exon-Intron Border Signal Scan')
plt.legend()
plt.grid(True, alpha=0.5)
plt.show()