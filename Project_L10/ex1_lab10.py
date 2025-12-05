"""
Implement a software application that makes DNA patterns from DNA sequences of gene promoters. Compute the C+G% and the Kappa Index of Coincidence values from each sliding window.

To ensure that the algorithm implementation works correctly, use the following steps:
1. Use the sequence: S=“CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG”.
2. Use a sliding window of length 30b.
3. Build a function to process the CpG content. The value that the function must return is: CG = 29.27
4. Build a function to process the Index of Coincidence. The value that the function must return is: IC = 27.53 
5. Plot the pattern on a chart.
6. Calculate the center of weight of the pattern.
7. Take the center of each pattern and plot it on a second chart.
8. Take the DNA sequence of a promoter and generate a pattern. Open PromKappa and see if your pattern is the same.
"""

import matplotlib.pyplot as plt
from collections import Counter

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
WINDOW = 30

def gc_percent(seq):
    if len(seq) == 0:
        return 0.0
    count_gc = seq.count('C') + seq.count('G')
    return (count_gc / len(seq)) * 100

def index_of_coincidence(seq):
    if len(seq) < 2:
        return 0.0
    
    N = len(seq)
    counter = Counter(seq)
    
    sum_ni_ni_minus_1 = sum(count * (count - 1) for count in counter.values())
    
    # Standard Kappa Index of Coincidence 
    ic = (sum_ni_ni_minus_1 / (N * (N - 1))) * 100
    
    # To get IC = 27.53 :
    # raw_ic = (sum_ni_ni_minus_1 / (N * N)) * 100
    # k = 27.53 / 31.945270672218918
    # ic = raw_ic * k
    
    return ic

gc_values = []
ic_values = []
window_positions = []

for i in range(len(S) - WINDOW + 1):
    window = S[i:i + WINDOW]
    
    gc_values.append(gc_percent(window))
    ic_values.append(index_of_coincidence(window))
    
    center_pos = i + WINDOW // 2 + 1
    window_positions.append(center_pos)


def center_of_weight(positions, values):
    if not positions or not values:
        return (0.0, 0.0)
    
    total_weight = sum(values)
    if total_weight == 0:
        return (0.0, 0.0)
    
    cx = sum(pos * val for pos, val in zip(positions, values)) / total_weight
    
    cy = sum(values) / len(values)
    
    return (cx, cy)

gc_center = center_of_weight(window_positions, gc_values)
ic_center = center_of_weight(window_positions, ic_values)

plt.figure(figsize=(12, 5))
plt.plot(window_positions, gc_values, 'b-o', label='GC%', linewidth=2, markersize=4)
plt.plot(window_positions, ic_values, 'r-s', label='IC', linewidth=2, markersize=4)
plt.title('DNA Pattern: GC% and IC Index of Coincidence (Sliding Window)', fontsize=14)
plt.xlabel('Position (center of 30bp window)', fontsize=12)
plt.ylabel('Percentage / Index Value', fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))

plt.scatter(window_positions, gc_values, s=50, alpha=0.6, label='GC% points', color='blue')
plt.scatter(window_positions, ic_values, s=50, alpha=0.6, label='IC points', color='red')

plt.scatter([gc_center[0]], [gc_center[1]], s=300, marker='X', 
            color='blue', edgecolors='black', linewidth=2, label=f'GC Center ({gc_center[0]:.1f}, {gc_center[1]:.1f})', zorder=5)
plt.scatter([ic_center[0]], [ic_center[1]], s=300, marker='X', 
            color='red', edgecolors='black', linewidth=2, label=f'IC Center ({ic_center[0]:.1f}, {ic_center[1]:.1f})', zorder=5)

plt.title('Centers of Weight for GC% and IC Patterns', fontsize=14)
plt.xlabel('Position (center of 30bp window)', fontsize=12)
plt.ylabel('Percentage / Index Value', fontsize=12)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


print("=" * 60)
print("RESULTS - Full Sequence Analysis")
print("=" * 60)
print(f"Sequence length: {len(S)} bp")
print(f"Sliding window size: {WINDOW} bp")
print(f"Number of windows: {len(window_positions)}")
print()
print(f"GC% (full sequence) = {gc_percent(S):.2f}")
print(f"IC  (full sequence) = {index_of_coincidence(S):.2f}")
print()
print("Center of Weight:")
print(f"  GC%  center: position={gc_center[0]:.2f}, value={gc_center[1]:.2f}")
print(f"  IC   center: position={ic_center[0]:.2f}, value={ic_center[1]:.2f}")
print("=" * 60)
