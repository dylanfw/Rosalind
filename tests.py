import io
import unittest
import rosalind


class TestRosalind(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None

    def test_count_bases(self):
        dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
        expected_counts = "20 12 17 21"
        self.assertEqual(rosalind.count_bases(dna), expected_counts)

    def test_dna_to_rna(self):
        dna = "GATGGAACTTGACTACGTAAATT"
        expected_rna = "GAUGGAACUUGACUACGUAAAUU"
        self.assertEqual(rosalind.dna_to_rna(dna), expected_rna)

    def test_reverse_complement(self):
        dna = "AAAACCCGGT"
        expected_reverse_complement = "ACCGGGTTTT"
        self.assertEqual(rosalind.reverse_complement(dna), expected_reverse_complement)

    def test_hamming_distance(self):
        dna1 = "GAGCCTACTAACGGGAT"
        dna2 = "CATCGTAATGACGGCCT"
        expected_distance = 7
        self.assertEqual(rosalind.hamming_distance(dna1, dna2), expected_distance)

    def test_rna_to_amino_acids(self):
        rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        expected_sequence = "MAMAPRTEINSTRING"
        self.assertEqual(rosalind.rna_to_amino_acids(rna), expected_sequence)

    def test_parse_fasta(self):
        dataset = """>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT"""
        expected_data = {
            "Rosalind_6404": "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
            "Rosalind_5959": "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
            "Rosalind_0808": "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT",
        }
        self.assertEqual(rosalind._parse_fasta(dataset), expected_data)

    def test_greatest_gc_content(self):
        dataset = """>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT"""
        expected_output = "Rosalind_0808, 60.919540"
        self.assertEqual(
            rosalind.greatest_gc_content(io.StringIO(dataset)), expected_output
        )


if __name__ == "__main__":
    unittest.main()
