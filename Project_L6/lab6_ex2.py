# Choose 5 restriction enzymes that are used to digest (cut) one input sequnce. 
# Your simulation of the e.g. must show the resulting DNA segments on one column
import os
import math
import matplotlib.pyplot as plt

# ==== Paths ====
SEQ_PATH = '/Users/axdenis03/Documents/University /BioInformatics/Project_L6/ex2.fasta'

# ==== Functions ====
def read_fasta(file_path):
    """Read a FASTA file and return the DNA sequence (ACGT only)."""
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return "".join(c for c in sequence.upper() if c in "ACGT")

def revcomp(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def find_cut_positions(dna, site):
    """Find all cut positions of a site and its reverse complement in the DNA."""
    cut_positions = set()
    for s in [site, revcomp(site)]:
        start = 0
        while True:
            pos = dna.find(s, start)
            if pos == -1:
                break
            cut_positions.add(pos)
            start = pos + 1
    return sorted(cut_positions)

def digest_sequence(dna, site):
    """Return fragment lengths after digesting dna with a restriction site."""
    cut_positions = find_cut_positions(dna, site)
    if not cut_positions:
        return [len(dna)]
    fragments = []
    prev = 0
    for cut in cut_positions:
        fragments.append(cut - prev)
        prev = cut
    fragments.append(len(dna) - prev)
    return fragments

def size_to_position_bp_dynamic(bp, max_bp):
    """Convert fragment size to relative gel position dynamically."""
    return 1.0 - (math.log10(bp) / math.log10(max_bp))

def plot_gel(all_fragments):
    """Plot simulated gel electrophoresis for all fragments with dynamic scaling."""
    # Determine maximum fragment length to scale dynamically
    max_frag = max(length for _, lengths in all_fragments for length in lengths)
    max_frag = max(max_frag, 1000)  # minimum scaling threshold

    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_facecolor("black")

    # Ladder lane
    ladder_x, ladder_w = 0.05, 0.18
    ax.add_patch(plt.Rectangle((ladder_x, 0.05), ladder_w, 0.90, facecolor="#1f1f2e", edgecolor="#111111"))

    ladder_bps = [3000, 2000, 1500, 1200, 1000, 800, 600, 400, 200]
    max_ladder = max(max_frag, max(ladder_bps))
    for bp in ladder_bps:
        if bp > max_ladder:  # skip ladder bands larger than maximum fragment
            continue
        y = 0.05 + 0.90*size_to_position_bp_dynamic(bp, max_ladder)
        th = 0.007
        ax.add_patch(plt.Rectangle((ladder_x+0.01, y-th/2), ladder_w-0.02, th, facecolor="#e8e8f0"))
        if bp in [3000, 1500, 500]:
            ax.text(ladder_x-0.03, y, f"{bp} bp", va="center", ha="right", fontsize=10, fontweight="bold")
    ax.text(ladder_x+ladder_w/2, 0.965, "Ladder", ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Enzyme lanes
    lane_w, lane_gap = 0.18, 0.12
    lane_x = ladder_x + ladder_w + lane_gap
    for enzyme, frag_lengths in all_fragments:
        # Lane background
        ax.add_patch(plt.Rectangle((lane_x, 0.05), lane_w, 0.90, facecolor="#1f1f2e", edgecolor="#111111"))
        # Fragments
        for length in sorted(frag_lengths, reverse=True):
            y = 0.05 + 0.90*size_to_position_bp_dynamic(length, max_ladder)
            th = max(0.004, 0.016 - 0.000003*length)
            ax.add_patch(plt.Rectangle((lane_x+0.01, y-th/2), lane_w-0.02, th, facecolor="#e8e8f0", edgecolor="#cfcfe8"))
            ax.text(lane_x+lane_w+0.02, y, f"{length} bp", va="center", ha="left", fontsize=8, color="#000000")
        # Lane label
        ax.text(lane_x+lane_w/2, 0.965, enzyme, ha="center", va="bottom", fontsize=10, fontweight="bold")
        lane_x += lane_w + lane_gap

    plt.title("Simulated Restriction Digest Gel", color="white", fontsize=14)
    plt.tight_layout()
    plt.show()


# ==== Main ====
sequence = read_fasta(SEQ_PATH)
print("Original sequence length =", len(sequence), "bp")

dna = sequence  # use full sequence

# Define enzymes
enzymes = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "PstI": "CTGCAG"
}

# Digest DNA
all_fragments = []
for enzyme, site in enzymes.items():
    frag_lengths = digest_sequence(dna, site)
    print(f"{enzyme}: fragment lengths (bp): {frag_lengths}")
    all_fragments.append((enzyme, frag_lengths))

# Plot gel
plot_gel(all_fragments)