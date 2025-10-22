from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt

GENETIC_CODE = {
    "UUU": "Phe", "UUC": "Phe",
    "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",
    "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr",
    "CAU": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys",
    "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "UAA": "Stop", "UAG": "Stop", "UGA": "Stop"
}

def read_fasta_sequence(filepath: str) -> str:
    record = next(SeqIO.parse(filepath, "fasta"))
    sequence = str(record.seq).upper().replace("T", "U")
    return sequence


def get_codon_frequencies(sequence: str) -> Counter:
    codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
    codons = [c for c in codons if len(c) == 3]
    return Counter(codons)


def plot_top_codons(codon_counter: Counter, title: str):
    top_codons = codon_counter.most_common(10)
    codons, counts = zip(*top_codons)
    plt.figure(figsize=(8, 4))
    plt.bar(codons, counts, color="mediumseagreen")
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()


def top_amino_acids(codon_counter: Counter) -> list:
    amino_counts = Counter()
    for codon, count in codon_counter.items():
        aa = GENETIC_CODE.get(codon)
        if aa and aa != "Stop":
            amino_counts[aa] += count
    return amino_counts.most_common(3)

if __name__ == "__main__":
    covid_seq = read_fasta_sequence("covid.fasta")
    influenza_seq = read_fasta_sequence("influenza.fasta")

    covid_codons = get_codon_frequencies(covid_seq)
    influenza_codons = get_codon_frequencies(influenza_seq)

    plot_top_codons(covid_codons, "Top 10 Codons - COVID-19")
    plot_top_codons(influenza_codons, "Top 10 Codons - Influenza")

    covid_top = {c for c, _ in covid_codons.most_common(10)}
    flu_top = {c for c, _ in influenza_codons.most_common(10)}
    common = covid_top.intersection(flu_top)
    print("Common top codons between COVID-19 and Influenza:", ", ".join(common) or "None")
    print()

    print("Top 3 Amino Acids for COVID-19:")
    for aa, count in top_amino_acids(covid_codons):
        print(f" - {aa}: {count}")
    print()

    print("Top 3 Amino Acids for Influenza:")
    for aa, count in top_amino_acids(influenza_codons):
        print(f" - {aa}: {count}")
    print()
