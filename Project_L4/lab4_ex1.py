# RNA table ( I did it with AI for a better format )
GENETIC_CODE = {
    "UUU": "Phe", "UUC": "Phe",
    # Leucine
    "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    # Isoleucine
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",
    # Methionine (Start codon)
    "AUG": "Met",
    # Valine
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    # Serine
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser",
    # Proline
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    # Threonine
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    # Alanine
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    # Tyrosine
    "UAU": "Tyr", "UAC": "Tyr",
    # Histidine
    "CAU": "His", "CAC": "His",
    # Glutamine
    "CAA": "Gln", "CAG": "Gln",
    # Asparagine
    "AAU": "Asn", "AAC": "Asn",
    # Lysine
    "AAA": "Lys", "AAG": "Lys",
    # Aspartic acid
    "GAU": "Asp", "GAC": "Asp",
    # Glutamic acid
    "GAA": "Glu", "GAG": "Glu",
    # Cysteine
    "UGU": "Cys", "UGC": "Cys",
    # Tryptophan
    "UGG": "Trp",
    # Arginine
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    # Glycine
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    # Stop codons
    "UAA": "Stop", "UAG": "Stop", "UGA": "Stop"
}

def translate_rna_to_protein(rna_sequence: str) -> str:
    protein = []
    rna_sequence = rna_sequence.upper().replace("T", "U")  # we convert DNA to RNA if we need to 
    
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) < 3:
            break  
        amino_acid = GENETIC_CODE.get(codon, "")
        if amino_acid == "Stop":
            break
        if amino_acid:
            protein.append(amino_acid)
    
    return "-".join(protein)


if __name__ == "__main__":
    example_seq = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA" # example RNA 
    protein = translate_rna_to_protein(example_seq)
    print("RNA sequence:", example_seq)
    print("Protein sequence:", protein)
    print()

    dna_seq = "ATGCCTTAGGCT" # changes T to U   
    print("DNA sequence:", dna_seq)
    print("Protein sequence:", translate_rna_to_protein(dna_seq))
    print()

    mixed_case = "guuGCTtaaGCT" # it translates guu to GUU and it stops at taa(UAA)
    print("Mixed case sequence:", mixed_case)
    print("Protein sequence:", translate_rna_to_protein(mixed_case))
    print()

    stop_seq = "ACUUAA" # UAA means stop
    print("Sequence with early stop:", stop_seq)
    print("Protein sequence:", translate_rna_to_protein(stop_seq))
    print()

    invalid_seq = "AUGNNN"  # NNN is not a valid 
    print("Invalid sequence:", invalid_seq)
    print("Protein sequence:", translate_rna_to_protein(invalid_seq))
    print()
