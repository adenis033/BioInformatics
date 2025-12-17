# Download 10 influenza genomes adapt your app from the previous assignment
# in order to scan each genome for possible motives. For each
# genome make a chart that shows the signal with most likely locations original
# fundamental motives

import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

MOTIF = "AGCAAAAGCAGG"
BASE_PATH = "/Users/axdenis03/Documents/University /BioInformatics/Project_L12/exercise2/"

def load_local_genomes():
    records = []
    files = ["sequence.fasta"] + [f"sequence-{i}.fasta" for i in range(2, 11)]
    
    for filename in files:
        full_path = os.path.join(BASE_PATH, filename)
        if os.path.exists(full_path):
            for record in SeqIO.parse(full_path, "fasta"):
                records.append(record)
                break 
    return records

def calc_signal(seq, motif):
    s_len, m_len = len(seq), len(motif)
    signal = np.zeros(s_len - m_len + 1)
    for i in range(len(signal)):
        sub = seq[i : i + m_len]
        signal[i] = sum(1 for a, b in zip(sub, motif) if a == b) / m_len
    return signal

def plot_results(genomes):
    if not genomes: return
    
    fig, axes = plt.subplots(len(genomes), 1, figsize=(10, 2 * len(genomes)), constrained_layout=True)
    if len(genomes) == 1: axes = [axes]

    for i, rec in enumerate(genomes):
        seq = str(rec.seq).upper()
        y = calc_signal(seq, MOTIF)
        peaks = np.where(y >= 0.8)[0]

        axes[i].plot(y, color='#333333', linewidth=1)
        axes[i].plot(peaks, y[peaks], 'ro', markersize=2)
        axes[i].set_title(rec.id, loc='left', fontsize=9)
        axes[i].set_ylim(0, 1.1)
        axes[i].grid(alpha=0.3)
        if i == len(genomes) - 1: axes[i].set_xlabel("Position")

    plt.show()

if __name__ == "__main__":
    data = load_local_genomes()
    plot_results(data)