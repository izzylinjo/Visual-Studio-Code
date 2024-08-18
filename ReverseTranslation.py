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

protein_sequence = "MATWYLGLPTWYLWCCKLWYTIILTWYGYPWLPHKGYPWFDAAVNGYPWRQQWTYIIL**"

def find_suitible_codon(dna_sequence, protein_sequence, dna_sequence_position):
    amino_acid_position = int(dna_sequence_position/3)
    amino_acid = protein_sequence[amino_acid_position]
    specific_codons = ecoli_codon_dict[amino_acid]
    used_codon = dna_sequence[dna_sequence_position: dna_sequence_position + 3]
    codon_rarity = specific_codons.index(used_codon)

    if len(specific_codons) > codon_rarity + 1: 
         return specific_codons[codon_rarity + 1]
    return specific_codons[codon_rarity]

def construct_dna_fragment(dna_sequence, protein_sequence, dna_position): 
    dna_fragment = ""
    for i in range(0, 7, 3):
        dna_fragment += find_suitible_codon(dna_sequence, protein_sequence, dna_position + i)

    return dna_fragment
    

def edit_dna_sequence(dna_sequence, protein_sequence, dna_position):
    strand_one = dna_sequence[:dna_position]
    middle_strand = construct_dna_fragment(dna_sequence, protein_sequence, dna_position)
    strand_two = dna_sequence[dna_position+9:]
    return strand_one + middle_strand + strand_two

def reverseTranslate(protein_seq, MOST_COMMON_CODON=0):
    dna_seq = ""
    for amino_acid in protein_seq:
        dna_seq += ecoli_codon_dict[amino_acid][MOST_COMMON_CODON]
    return dna_seq

#copiolot addition need to test further but is a sliding window algorithm 
def get_repeating_fragments(dna_seq, fragmentSize=9):
    fragments = []
    repeats = {}
    for i in range(len(dna_seq) - fragmentSize+ 1):
        fragment = dna_seq[i:i + fragmentSize]
        if fragment in fragments:
            repeats[fragment] = i 
        fragments.append(fragment)
    return repeats

def adjust_to_multiple_of_threee(postion):
    if postion%3 == 1:
        return postion - 1
    elif postion%3 == 2:
        return postion + 1
    return postion

def position_of_repeats(fragments):
    dna_positions = [adjust_to_multiple_of_threee(i) for i in fragments.values()]
    sorted_dna_positions = []
    [sorted_dna_positions.append(x) for x in dna_positions if x not in sorted_dna_positions]
    return sorted_dna_positions

def loop_through_edit(dna_sequence, protein_sequence, ordered_positions):
    for position in ordered_positions:
        dna_sequence = edit_dna_sequence(dna_sequence, protein_sequence, position)
    return dna_sequence

def count_codon_differences(seq1, seq2):
    # Ensure both sequences are of the same length by padding the shorter one with spaces
    max_len = max(len(seq1), len(seq2))
    seq1 = seq1.ljust(max_len)
    seq2 = seq2.ljust(max_len)
    
    # Count differing codons
    differences = 0
    for i in range(0, max_len, 3):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]
        if codon1 != codon2:
            differences += 1
    return differences



def check_gc_content():
    pass

dna_sequence = reverseTranslate(protein_sequence)
fragments = get_repeating_fragments(dna_sequence)
ordered_positions = position_of_repeats(fragments)
edited_dna = loop_through_edit(dna_sequence, protein_sequence, ordered_positions)

codon_differences = count_codon_differences(dna_sequence, edited_dna)
print(f"Codon Differences: {codon_differences}")

print(edited_dna)
