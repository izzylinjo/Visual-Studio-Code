import unittest

class TestReverseTranslation(unittest.TestCase):
    def test_gc_content(self):
        dna_seq = aa_to_dna(aa_seq, ecoli_codon_dict)
        for i in range(len(dna_seq) - 11):
            self.assertTrue(30 < gc_content(dna_seq[i:i+12]) < 70, "GC content rule violated")

    def test_no_repeats(self):
        dna_seq = aa_to_dna(aa_seq, ecoli_codon_dict)
        self.assertFalse(has_repeats(dna_seq), "No repeat rule violated")

    def test_codon_usage(self):
        dna_seq = aa_to_dna(aa_seq, ecoli_codon_dict)
        for aa, codons in ecoli_codon_dict.items():
            used_codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3) if dna_seq[i:i+3] in codons]
            self.assertTrue(used_codons[0] == codons[0], "Codon usage rule violated")

if __name__ == '__main__':
    unittest.main()