# Find in seq S only the din and tri that exists without the use of the BF engine .
# In order to achieve the results one must verify this comb starting from beggining
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

dinucleotides = {}
for i in range(len(S) - 1):
    chunk = S[i:i+2]
    if chunk not in dinucleotides:
        dinucleotides[chunk] = 1

trinucleotides = {}
for i in range(len(S) - 2):
    chunk = S[i:i+3]
    if chunk not in trinucleotides:
        trinucleotides[chunk] = 1

print(dinucleotides)
print(trinucleotides)