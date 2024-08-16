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

def translate(dnaSeq): 
    pass 

def reverseTranslate(proteinSeq, i=0):
    dnaSeq = ""
    for x in proteinSeq:
        dnaSeq += ecoli_codon_dict[x][i]
    return dnaSeq
def checkGCContent():
    pass

def getRepition(dnaSeq):
    print(dnaSeq[:9])

dnaSeq = reverseTranslate(proteinSequence)

getRepition(dnaSeq)

