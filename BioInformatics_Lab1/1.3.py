def analyze_fasta(file_path):

    sequence = []

    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line.startswith(">"):  # skip FASTA header lines
                    sequence.append(line)
    except FileNotFoundError:
        print(f"ERROR: File '{file_path}' not found.")
        return

    seq = "".join(sequence).upper()

    if not seq:
        print(f"ERROR: No sequence found in file '{file_path}'.")
        return

    alphabet = sorted(set(seq))

    print(f"Analysis for file: {file_path}")
    print(f"Sequence length: {len(seq)} characters")
    print("-------------------------------------------\n")
    print(f"Alphabet of the sequence ({len(alphabet)} unique characters):")
    print(" ".join(alphabet))
    print("\nRelative character frequency:")

    total = len(seq)
    for ch in alphabet:
        count = seq.count(ch)
        percent = (count / total) * 100
        print(f"{ch} -> {count} occurrences ({percent:.2f}%)")

file_name = "sequence.fasta"
analyze_fasta(file_name)
