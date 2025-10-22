# Desing an app that uses the sliding window method
# in order to read the tm over the seq S.
# Use a sliding window f 8 positions and choose a FASTA file for input

import math

filename = input("Enter the full path of your FASTA file: ").strip()
if not filename:
    print("No file provided. Exiting.")
    exit()

sequence = ""
with open(filename, 'r') as f:
    for line in f:
        if not line.startswith(">"):
            sequence += line.strip().upper()

window_size = 8
Na = 0.0001

print(f"{'Window':>6} {'Seq':>10} {'GC%':>6} {'Tm1':>6} {'Tm2':>6}")
print("-" * 40)

for i in range(len(sequence) - window_size + 1):
    window = sequence[i:i+window_size]

    A = window.count('A')
    T = window.count('T')
    G = window.count('G')
    C = window.count('C')

    gc = ((G + C) / window_size) * 100

    tm1 = 4 * (G + C) + 2 * (A + T)
    tm2 = 81.5 + 16.6 * math.log10(Na) + 0.41 * gc - (600 / window_size)

    print(f"{i+1}-{i+window_size:<3} {window:>10} {gc:6.2f} {tm1:6.2f} {tm2:6.2f}")
