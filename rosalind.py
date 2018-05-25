"""
A collection of various methods for solving bioinformatics challenges from Rosalind.info
"""


def count_bases(dna):
    """Return the number of Adenine, Cytosine, Guanine, and Thymine bases in a given DNA string"""
    a_count, c_count, g_count, t_count = 0, 0, 0, 0
    for base in dna.upper():
        if base == "A":
            a_count += 1
        elif base == "C":
            c_count += 1
        elif base == "G":
            g_count += 1
        elif base == "T":
            t_count += 1
        else:
            raise ValueError(f"Unknown nucleotide found: {base}")
    return f"{a_count} {c_count} {g_count} {t_count}"


def dna_to_rna(dna):
    """Translate the given DNA string into RNA"""
    return dna.translate(str.maketrans("T", "U"))


def reverse_complement(dna):
    """Return the reverse complement of the given DNA string"""
    reversed_dna = dna[::-1]
    return reversed_dna.translate(str.maketrans("ACGT", "TGCA"))


def greatest_gc_content(fasta_file):
    """Given a file in FASTA format, return the ID and GC content of the DNA string with the highest GC content"""
    data = _parse_fasta(fasta_file.read())
    gc_data = [
        {"name": name, "string": string, "gc": _gc_content(string)}
        for name, string in data.items()
    ]
    greatest_gc = max(gc_data, key=lambda d: d["gc"])
    return "{}, {:0.6f}".format(greatest_gc["name"], greatest_gc["gc"])


def _gc_content(genetic_string):
    """Return the percentage of G and C nucleotides in the given genetic string"""
    gc_count = genetic_string.count("G") + genetic_string.count("C")
    return gc_count / len(genetic_string) * 100


def hamming_distance(dna1, dna2):
    """Return the number of differing bases between two strings of DNA"""
    return sum([1 for i, base in enumerate(dna1) if base != dna2[i]])


def rna_to_amino_acids(rna):
    """Translate the given RNA string into a sequence of amino acids using the RNA_CODON_TABLE"""
    rna = rna.upper()
    amino_acids = []
    for i in range(0, len(rna), 3):
        codon = rna[i : i + 3]
        amino_acid = RNA_CODON_TABLE.get(codon)
        if not amino_acid:
            continue
        if amino_acid == "*":
            break
        amino_acids.append(amino_acid)
    return "".join(amino_acids)


def _parse_fasta(fasta):
    data = {}
    records = fasta.split(">")

    # First line should begin with >, so we eliminate the empty element
    records = records[1:]
    for record in records:
        record = record.split("\n", 1)
        try:
            name = record[0].strip()
            genetic_string = "".join(record[1].split("\n"))
            data[name] = genetic_string
        except KeyError:
            raise ValueError("Unable to parse FASTA entry")
    return data


RNA_CODON_TABLE = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
