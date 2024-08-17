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

proteinSequence = "MATWYLGLPTWYLWCCKLWYTIILTWYGYPWLPHKGYPWFDAAVNGYPWRQQWTYIIL**"

def translate(dna_seq): 
    pass 

def find_suitible_codon():

    pass 

def reverseTranslate(protein_seq, MOST_COMMON_CODON=0):
    dnaSeq = ""
    for amino_acid in protein_seq:
        dnaSeq += ecoli_codon_dict[amino_acid][MOST_COMMON_CODON]
    return dnaSeq

#copiolot addition need to test further but is a sliding window algorithm 
def get_repeating_fragments(dna_seq, fragmentSize=9):
    repeating_fragments = []
    fragments = []
    for i in range(len(dna_seq) - fragmentSize+ 1):
        fragment = dna_seq[i:i + fragmentSize]
        if fragment in fragments:
            repeating_fragments.append(fragment)
        fragments.append(fragment)
    return repeating_fragments

def check_gc_content():
    pass

def getRepition(dnaSeq):

    pass

dna_sequence = reverseTranslate(proteinSequence)
print(dna_sequence)

fragments = get_repeating_fragments(dna_sequence)
print(fragments)