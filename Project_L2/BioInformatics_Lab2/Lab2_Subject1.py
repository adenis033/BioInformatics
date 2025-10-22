# Find the percentage for all the dinucleotide and trinucleotide combinations for the sequence:
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

# 1. Build a brute force engine to generate all dinucleotide and trinucleotide combinations.

alph = sorted(set(S))

dinucleotides = []

for first in alph:
    for second in alph:
        dinucleotides.append(first + second)

print(dinucleotides)

trinucleotide = []

for first in alph:
    for second in alph:
        for third in alph:
            trinucleotide.append(first + second + third)

print(trinucleotide)

print()
# 2. For each combination, find out the percentage inside the S sequence.

countDin = {}
for i in range(len(S)-1):
    chunk = S[i:i+2]
    if chunk in countDin:
        countDin[chunk] += 1
    else:
        countDin[chunk] = 1

print(countDin)

countTri = {}
for x in range(len(S)-2):
    chunk = S[x:x+3]
    if chunk in countTri:
        countTri[chunk] += 1
    else:
        countTri[chunk] = 1

print(countTri)

# 3. Show the percentage for each combination in the output of your implementation.

print()

dinTotal = sum(countDin.values())
for chunk in countDin:
    print(chunk,"-",round(countDin[chunk]/dinTotal*100, 2), "%")

print()

triTotal = sum(countTri.values())
for chunk in countTri:
    print(chunk,"-",round(countTri[chunk]/triTotal*100, 2), "%")

