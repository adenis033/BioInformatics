import random

def random_dna(n):
    return "".join(random.choice("ACGT") for _ in range(n))

def insert_transposon(seq, transposon, pos):
    return seq[:pos] + transposon + seq[pos:]

def find_transposons(seq, patterns):
    results = []
    for p in patterns:
        start = 0
        l = len(p)
        while True:
            i = seq.find(p, start)
            if i == -1:
                break
            results.append((p, i, i + l))
            start = i + 1
    return results

dna = random_dna(random.randint(200, 400))

transposons = [
    "ATGCGTATGCGT",
    "GGATCCGGATCC",
    "TTACGTTTACGT",
    "CAGTCAGT"
]

for t in transposons:
    pos = random.randint(0, len(dna))
    dna = insert_transposon(dna, t, pos)

hits = find_transposons(dna, transposons)

print("Final DNA:\n", dna)
print()
print("Detected transposons:")
for h in hits:
    print(h)