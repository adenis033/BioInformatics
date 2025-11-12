import random
import os

fasta_path = '/Users/axdenis03/Documents/University /BioInformatics/Project_L5/sequence.fasta'

def read_fasta(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith(">"))
    return sequence.upper()

dna_sequence = read_fasta(fasta_path)

# a & b
samples = []
for _ in range(2000):
    length = random.randint(100, 150)
    start = random.randint(0, len(dna_sequence) - length)
    sample = dna_sequence[start:start + length]
    samples.append(sample)

# c
min_overlap = 10
reads = samples.copy()
current_seq = reads.pop(0)

merged = True
while merged and reads:
    merged = False
    best_overlap = 0
    best_index = None
    attach_to_back = True  # True = back, False = front

    for i, read in enumerate(reads):
        max_overlap = min(len(current_seq), len(read))

        # check overlap at the back
        for o in range(min_overlap, max_overlap + 1):
            if current_seq[-o:] == read[:o] and o > best_overlap:
                best_overlap = o
                best_index = i
                attach_to_back = True

        # check overlap at the front
        for o in range(min_overlap, max_overlap + 1):
            if read[-o:] == current_seq[:o] and o > best_overlap:
                best_overlap = o
                best_index = i
                attach_to_back = False

    if best_index is not None:
        read = reads.pop(best_index)
        if attach_to_back:
            current_seq += read[best_overlap:]
        else:
            current_seq = read[:-best_overlap] + current_seq
        merged = True

# truncate to original length to avoid going over
reconstructed_sequence = current_seq[:len(dna_sequence)]

# d
output_dir = os.path.dirname(fasta_path)

recon_path = os.path.join(output_dir, "reconstructed_sequence.txt")
with open(recon_path, "w") as f:
    f.write(reconstructed_sequence)

dna_path = os.path.join(output_dir, "dna_sequence.txt")
with open(dna_path, "w") as f:
    f.write(dna_sequence)

answer_path = os.path.join(output_dir, "answer.txt")
explanation = """
The algorithm tries to rebuild the DNA sequence by merging random short reads when they overlap for at least 10 basis/letters.

problems with this approach:
	1.	Overlaps in the middle are missed:
	it only tries to merge reads at the current sequence’s start or end, ignoring overlaps in the middle.
	if a read overlaps somewhere in the middle, it will not merge.
	2.	Order matters:
	reads are merged in the order they appear.
	if multiple reads could fit in different places, the sequence may be reconstructed incorrectly.
	3.	Incomplete reconstruction:
	some regions may not be included in any read, so they cannot be reconstructed.
	the reconstructed sequence can end up shorter than the original.
    more reads more coverage, better chance of full reconstruction.
"""
with open(answer_path, "w") as f:
    f.write(explanation.strip())

print("dna_sequence.txt created")
print("reconstructed_sequence.txt created")
print("answer.txt created")