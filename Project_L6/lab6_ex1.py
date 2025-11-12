'''
Gel electrophoresis is an analysis method implemented in all disciplines of life sciences. 
The results of gel electrophoresis indicate the relative sizes of fragments, 
which is useful for restriction mapping and analyzing PCR fragments. 
1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
2. Take 10 random samples from this sequence, between 100-3000 bases.
3. Store these samples in an array.
4. Simulate the migration of these DNA segments on the electrophoresis gel, based on their molecular weights - however, 
their length should be sufficient for this exercise (show a visual representation).
Note: Short DNA fragments meet small friction forces and travel faster through the electrophoresis gel. 
Long DNA fragments exhibit a high friction force and travel slowly through the electrophoresis gel.
'''

import random
import matplotlib.pyplot as plt

def read_fasta(file_path):
    dna_sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                dna_sequence += line.strip()
    return dna_sequence.upper()

file_path = '/Users/axdenis03/Documents/University /BioInformatics/Project_L6/sequence.fasta'
original_sequence = read_fasta(file_path)
sequence_length = len(original_sequence)

num_samples = 10
fragments = []
for _ in range(num_samples):
    frag_length = random.randint(100, min(3000, sequence_length))
    start_pos = random.randint(0, sequence_length - frag_length)
    fragments.append(original_sequence[start_pos:start_pos + frag_length])

fragments_bp = [len(frag) for frag in fragments]
ladder_bp = [3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]

all_bps = fragments_bp + ladder_bp
max_bp, min_bp = max(all_bps), min(all_bps)
MAX_MIGRATION, MIN_MIGRATION = 10, 1

def migration(length):
    normalized = (length - min_bp) / (max_bp - min_bp)
    return MIN_MIGRATION + (1 - normalized) * (MAX_MIGRATION - MIN_MIGRATION)

sample_distances = [migration(l) for l in fragments_bp]
ladder_distances = [migration(l) for l in ladder_bp]

plt.figure(figsize=(4, 7))
ax = plt.gca()
ax.set_facecolor('black')

ax.hlines(y=ladder_distances, xmin=0.65, xmax=1.35, color='white', linewidth=3)
ax.hlines(y=sample_distances, xmin=1.65, xmax=2.35, color='white', linewidth=3)

for bp, dist in zip(ladder_bp, ladder_distances):
    if bp in [3000, 1500, 500]:
        plt.text(0.55, dist, f"{bp} bp -", va='center', ha='right', fontsize=10, color='white')

plt.xlim(0, 3)
plt.ylim(0, MAX_MIGRATION + 1)
ax.invert_yaxis()
plt.xticks([1, 2], ["Ladder", "Samples"], color='white')
plt.yticks([])
plt.title("Simulated Gel Electrophoresis", color='white')
for side in ['top', 'right', 'left']:
    ax.spines[side].set_visible(False)
ax.spines['bottom'].set_color('white')
ax.tick_params(axis='x', length=0, colors='white')
plt.tight_layout()
plt.show()