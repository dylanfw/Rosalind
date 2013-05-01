from string import maketrans

def dnaBases(dna):
    # counts the occurences of A, C, G, and T in a given DNA strand
    adenine = 0
    cytosine = 0
    guanine = 0
    thymine = 0
    for base in dna.upper():
        if base == "A":
            adenine += 1
        elif base == "C":
            cytosine += 1
        elif base == "G":
            guanine += 1
        elif base == "T":
            thymine += 1
        else:
            return "Unknown nucleotide given: " + base
    return "%d %d %d %d" % (adenine, cytosine, guanine, thymine)

def dnaToRna(dna):
    # translate a given DNA strand into RNA
    trans = maketrans("T", "U")
    return dna.translate(trans)

def reverseComplement(dna):
    # returns the reverse complement to a given DNA strand
    reverseDNA = dna[::-1]
    trans = maketrans("ACGT", "TGCA")
    return reverseDNA.translate(trans)

def greatestGC(listFile):
    # NEEDS REWORKED!
    # given a file containing a list of proper names and dna strings, returns
    # the name and value of the dna string with the highest GC-percentage (#ofG + #ofC / totalbases)
    names = []
    gc = []
    total = 0
    with open(listFile) as f:
        data = f.read()
    entries = data.split(">")
    entries.remove("")
    for entry in entries:
        name = entry[:entry.find("\n")]
        names.append(name)
        dna = entry[entry.find(name)+len(name):]
        dna = dna.split("\n")
        dna = "".join(dna)
        entry_gc = 0.0
        for base in dna:
            if base == "G" or base == "C":
                entry_gc += 1
        gc.append((entry_gc/len(dna))*100)
    greatestVal = 0
    for x in range(0, len(gc)):
        if gc[x] > greatestVal:
            greatestVal = gc[x]
            greatestName = names[x]
    return greatestName, greatestVal

def hammingDistance(dna1, dna2):
    # compares two given dna strings, returns the number of differing bases
    dh = 0
    for x in range(0, len(dna1)):
        if dna1[x] != dna2[x]:
            dh += 1
    return dh

def replace_all(text, dic):
    # uses provided dictionary as a mapping to translate all occurences in text
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text

def proteinTranslation(rna):
    # given a string of rna, return the corresponding sequence of protein translations
    rna = rna.upper()
    acids = []
    table = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S',
            'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
            'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L',
            'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R',
            'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N',
            'AAA': 'K', 'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A',
            'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    i = 0
    while(i < len(rna)):
        for codon, acid in table.iteritems():
            if rna[i:i+3] == codon:
                acids.append(acid)
        i += 3
    return "".join(acids)