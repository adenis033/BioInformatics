import numpy as np
import matplotlib.pyplot as plt

fasta_path = '/Users/axdenis03/Documents/University /BioInformatics/Project_L9/sequence.fasta'

def read_fasta(path):
    seq = ""
    try:
        with open(path, 'r') as f:
            for line in f:
                if not line.startswith(">"):
                    seq += line.strip().upper()
    except FileNotFoundError:
        print(f"Error: File not found at {path}")
        return ""
    return seq

dna = read_fasta(fasta_path)

if not dna:
    print("No DNA sequence loaded. Check file path.")
    exit()

enzymes = {
    "EcoRI":  ("GAATTC", 1),
    "BamHI":  ("GGATCC", 1),
    "HindIII":("AAGCTT", 1),
    "TaqI":   ("TCGA",   1),
    "HaeIII": ("GGCC",   2)
}

def digest(seq, recog, cut_offset):
    sites = []
    start = 0
    while True:
        idx = seq.find(recog, start)
        if idx == -1:
            break
        sites.append(idx + cut_offset)
        start = idx + 1
    
    if not sites:
        return [], [len(seq)]
        
    sites.sort()
    fragments = []
    prev = 0
    for site in sites:
        fragments.append(site - prev)
        prev = site
    fragments.append(len(seq) - prev)
    
    return sites, sorted(fragments, reverse=True)

lane_results = {}
print(f"DNA Length: {len(dna)} bp\n")

for name, (seq_pat, offset) in enzymes.items():
    cuts, frags = digest(dna, seq_pat, offset)
    lane_results[name] = frags
    print(f"--- {name} ---")
    print(f"Number of cleavages: {len(cuts)}")
    print(f"Cleavage positions: {cuts}")
    print(f"Fragment lengths: {frags}\n")

marker_sizes = [3000, 2500, 2000, 1500, 1200, 1000, 800, 600, 400, 200]
lane_results["Marker"] = marker_sizes

gel_h = 1000
gel_w = len(lane_results) * 100
gel_img = np.zeros((gel_h, gel_w))

min_bp = 100
max_bp = 3500

def bp_to_y(bp, h_px, min_b, max_b):
    bp = max(min_b, min(bp, max_b))
    rel_pos = (np.log10(max_b) - np.log10(bp)) / (np.log10(max_b) - np.log10(min_b))
    return int(rel_pos * (h_px - 50)) + 25

col_width = 100
band_width = 60
keys = list(lane_results.keys())

for idx, name in enumerate(keys):
    fragments = lane_results[name]
    center_x = idx * col_width + (col_width // 2)
    
    for frag in fragments:
        y = bp_to_y(frag, gel_h, min_bp, max_bp)
        
        y_start = max(0, y-2)
        y_end = min(gel_h, y+3)
        x_start = center_x - (band_width // 2)
        x_end = center_x + (band_width // 2)
        gel_img[y_start:y_end, x_start:x_end] = 1.0 
        
        y_s_start = max(0, y-6)
        y_s_end = min(gel_h, y+7)
        gel_img[y_s_start:y_s_end, x_start:x_end] += 0.3

plt.figure(figsize=(10, 8), facecolor='black')
plt.imshow(gel_img, cmap='gray', aspect='auto', vmin=0, vmax=2)
plt.axis('off')

for idx, name in enumerate(keys):
    x = idx * col_width + (col_width // 2)
    plt.text(x, -20, name, color='white', ha='center', fontsize=11, rotation=45)

for size in marker_sizes:
    y = bp_to_y(size, gel_h, min_bp, max_bp)
    plt.text(gel_w + 10, y, f"{size}", color='white', va='center', fontsize=9)

plt.text(gel_w + 10, -20, "bp", color='white', ha='left', fontsize=11)
plt.title(f"Restriction Digest: {len(dna)} bp", color='white', pad=30)
plt.tight_layout()
plt.show()