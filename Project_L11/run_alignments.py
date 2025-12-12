# Download 10 HIV-1 genomes and use them as an input for your application
# Adapt your current algorithm for peer-wise alignment of 2 genomes.
# Note *take into consideration the fact that your current implementation
# can not align 2 genomes directly, therefore, you must find a strategy in which
# the 2 genomes can be aligned step-by-step on smaller regions

# Make the alignment between each genome

import os
from itertools import combinations

BASE_DIR = "/Users/axdenis03/Documents/University /BioInformatics/Project_L11"

MATCH = 1
MISMATCH = -1
GAP = -1
WINDOW = 800
OVERLAP = 150
SLACK = 300
STEP = 50

def read_fasta(path):
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
    return "".join(lines[1:]).upper()

def needleman_wunsch(s1, s2):
    n, m = len(s1), len(s2)
    M = [[0]*(m+1) for _ in range(n+1)]
    T = [[""]*(m+1) for _ in range(n+1)]

    for i in range(1, n+1):
        M[i][0] = i*GAP; T[i][0] = "U"
    for j in range(1, m+1):
        M[0][j] = j*GAP; T[0][j] = "L"

    for i in range(1, n+1):
        for j in range(1, m+1):
            d = M[i-1][j-1] + (MATCH if s1[i-1]==s2[j-1] else MISMATCH)
            u = M[i-1][j] + GAP
            l = M[i][j-1] + GAP
            best = max(d, u, l)
            M[i][j] = best
            T[i][j] = "D" if best==d else ("U" if best==u else "L")

    a1, a2, mid = "", "", ""
    i, j = n, m
    matches = 0
    while i>0 or j>0:
        if i>0 and j>0 and T[i][j]=="D":
            c1, c2 = s1[i-1], s2[j-1]
            a1 = c1 + a1; a2 = c2 + a2
            mid = ("|" if c1==c2 else ".") + mid
            if c1==c2: matches+=1
            i-=1; j-=1
        elif i>0 and T[i][j]=="U":
            a1 = s1[i-1] + a1; a2 = "-" + a2; mid = " " + mid; i-=1
        else:
            a1 = "-" + a1; a2 = s2[j-1] + a2; mid = " " + mid; j-=1

    similarity = round(matches/len(a1)*100,2)
    return a1, mid, a2, matches, similarity

def find_best_match(window, B, start, end):
    best_sim = -1
    best = None
    for pos in range(start, end, STEP):
        sub = B[pos:pos+len(window)+200]
        if not sub: continue
        a1, mid, a2, m, sim = needleman_wunsch(window, sub)
        if sim > best_sim:
            best_sim = sim
            best = (a1, mid, a2, pos)
    # fallback to global alignment if nothing found
    if best is None:
        best = needleman_wunsch(window, B[start:end]) + (start,)
    return best

def stitch(stitched, new_aln):
    if stitched is None:
        return new_aln[0], new_aln[1], new_aln[2]
    s1, mid1, s2 = stitched
    a1, mid, a2 = new_aln[:3]
    return s1 + a1, mid1 + mid, s2 + a2

def align_genomes(A, B, out_file):
    windows = [A[i:i+WINDOW] for i in range(0, len(A), WINDOW-OVERLAP)]
    stitched = None
    prev_end = 0

    for idx, w in enumerate(windows):
        start = max(0, prev_end - SLACK)
        end = min(len(B), prev_end + SLACK + len(w))
        aln = find_best_match(w, B, start, end)
        a1, mid, a2, pos = aln[:4]
        prev_end = pos + len(a2.replace("-", ""))
        stitched = stitch(stitched, (a1, mid, a2))

    if stitched is None:
        stitched = needleman_wunsch(A, B)

    top, midline, bottom = stitched
    matches = sum(1 for x,y in zip(top,bottom) if x==y)
    sim = round(matches/len(top)*100,2)

    with open(out_file,"w") as f:
        f.write("ALIGNMENT\n")
        f.write(top+"\n")
        f.write(midline+"\n")
        f.write(bottom+"\n\n")
        f.write(f"Matches: {matches}\n")
        f.write(f"Length: {len(top)}\n")
        f.write(f"Similarity: {sim}%\n")

def main():
    files = []
    for i in range(1,11):
        name = "sequence.fasta" if i==1 else f"sequence-{i}.fasta"
        full = os.path.join(BASE_DIR,name)
        if os.path.exists(full):
            files.append(full)

    if not files:
        print("No FASTA files found. Check your BASE_DIR!")
        return

    seqs = {f: read_fasta(f) for f in files}
    results_dir = os.path.join(BASE_DIR, "results")
    os.makedirs(results_dir, exist_ok=True)

    print("\nStarting pairwise alignments...\n")
    for f1, f2 in combinations(files,2):
        print("Aligning:", os.path.basename(f1), "vs", os.path.basename(f2))
        A = seqs[f1]; B = seqs[f2]
        out = os.path.join(results_dir, os.path.basename(f1)+"_vs_"+os.path.basename(f2)+".txt")
        align_genomes(A,B,out)

    print("\nAll alignments done! Check the 'results' folder.")

if __name__=="__main__":
    main()