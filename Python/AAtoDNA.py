# 1. few/no repeat nt longer than 8 bp 
# 2. 30 <  GC < 70 for all 12 nt sub seq 
# 3a. Always use most frequent codon 
# 3b. Be harmonious as in percentage

# No fragment reptition such as 'AAA AAA AAG AAA AAA AAG' 

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

p = "MIMAAALLLAAALLLAAALLLAKLAKLAKLAALLAALLCL*"

plist = list(p)
ecoliAA = list(ecoli_codon_dict.keys())
isProtein = all(x in ecoliAA for x in plist)
#convert amino acids to corresponding DNA Sequence regardless of repetition 
def aatoDNA(aa):
    dna = ""
    if isProtein:
        for x in aa:
            dna = dna + ecoli_codon_dict[x][0]
        return dna
    else: 
        return "Not Valid Amino Acid Sequence."
    

def checkAARep(aaSeq):
    
    return False     
   
# retroactively check if fragments are reptiting
def checkFragments(dnaSeq, fraglength = 9, move = 0): 
    fragments = [(dnaSeq[i: i+fraglength]) for i in range(move, len(dnaSeq), fraglength)] 

    for i in range(0, len(fragments)): 
       if fragments[i] == fragments[i-1]:
          x=2 

    if fraglength < len(dnaSeq)/2:
        return checkFragments(dnaSeq, fraglength + 3)

    
    # return editFragments(dnaSeq, list of posititons needed) 
    

# edit repetitive fragments  
def editFragment(dnaSeq, positions): 
    #positions and edit them out 
    return False

print(checkFragments(aatoDNA(plist)))


#def checkGC(): 

    
#def editGC(): 



