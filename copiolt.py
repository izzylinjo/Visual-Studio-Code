ecoli_codon_dict = {
    'A': ['GCG', 'GCC', 'GCA', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGC', 'GGT', 'GGG', 'GGA'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'K': ['AAA', 'AAG'],
    'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCG', 'CCA', 'CCT', 'CCC'],
    'Q': ['CAG', 'CAA'],
    'R': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'], 
    'S': ['AGC', 'TCT', 'AGT', 'TCC', 'TCA', 'TCG'],
    'T': ['ACC', 'ACG', 'ACT', 'ACA'],
    'V': ['GTG', 'GTT', 'GTC', 'GTA'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TGA', 'TAG']
}

p = "MATWYLGLPTWYLWCCKLWYTIILTWYGYPWLPHKGYPWFDAAVNGYPWRQQWTYIIL**"
#p2 = "MIMAAALAAALMMLAMMLAPPPLPPPLPPPPPPLLLLLL"

plist = list(p)

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

def has_repeats(seq, length=8):
    for i in range(len(seq) - length + 1):
        if seq[i:i+length] in seq[i+length:]:
            return True
    return False

def find_suitable_codon(dna_seq, codons):
    for codon in codons:
        temp_seq = dna_seq + codon
        if not has_repeats(temp_seq) and 30 < gc_content(temp_seq[-12:]) < 70:
            return codon
    return codons[-1]

def aa_to_dna(aa_seq, ecoli_codon_dict):
    dna_seq = ''
    for aa in aa_seq:
        codon = find_suitable_codon(dna_seq, ecoli_codon_dict[aa])
        dna_seq += codon
    return dna_seq

dna_seq = aa_to_dna(p, ecoli_codon_dict)
print(dna_seq)

import unittest

class TestReverseTranslation(unittest.TestCase):
    def test_gc_content(self):
        dna_seq = aa_to_dna(p, ecoli_codon_dict)
        for i in range(len(dna_seq) - 11):
            self.assertTrue(30 < gc_content(dna_seq[i:i+12]) < 70, "GC content rule violated")

    def test_no_repeats(self):
        dna_seq = aa_to_dna(p, ecoli_codon_dict)
        self.assertFalse(has_repeats(dna_seq), "No repeat rule violated")

    def test_codon_usage(self):
        dna_seq = aa_to_dna(p, ecoli_codon_dict)
        for aa, codons in ecoli_codon_dict.items():
            used_codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3) if dna_seq[i:i+3] in codons]
            self.assertTrue(used_codons[0] == codons[0], "Codon usage rule violated")

if __name__ == '__main__':
    unittest.main()