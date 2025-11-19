from Bio import SeqIO

def revcomp(s):
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join(comp[b] for b in reversed(s))

def find_inverted_repeats(seq, min_len=4, max_len=6):
    hits = []
    n = len(seq)
    for k in range(min_len, max_len+1):
        for i in range(n-k):
            left = seq[i:i+k]
            right = revcomp(left)
            j = seq.find(right, i+k)
            if j != -1:
                hits.append((left, i, j, j+k))
    return hits

genomes = ["genome1.fasta","genome2.fasta","genome3.fasta"]

for g in genomes:
    record = SeqIO.read(g, "fasta")
    seq = str(record.seq)
    results = find_inverted_repeats(seq)
    print("\nGenome:", g)
    print("Total inverted repeats found:", len(results))
    print("First 20 hits:")
    for r in results[:20]:
        print(r)