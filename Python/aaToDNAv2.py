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

p = "MIMAAALLLAAALLLAAALLLAKLAKLAKLAALLAALLL*"
#p2 = "MIMAAALAAALMMLAMMLAPPPLPPPLPPPPPPLLLLLL"

plist = list(p)
ecoliAA = list(ecoli_codon_dict.keys())

try:
    assert all(x in ecoliAA for x in plist)
except:
    exit('Seq does not contain all AAs')


#convert amino acids to corresponding DNA Sequence regardless of repetition 
def aatoDNA(aa, i=0):
    dna = ""
    for x in aa:
        dna += ecoli_codon_dict[x][i]
    return dna
    
def fragmentOccurence(dnaFragments, dnaSeq): 
        fragmentCount = []

        for x in dnaFragments: 
            fragmentCount.append(dnaSeq.count(x))
 
        return fragmentCount
      
    
def slidingFragmentation(dnaSeq, fraglength=3, shift=0 ):
    slidingFragments = [(dnaSeq[i: i+fraglength]) for i in range(shift, len(dnaSeq), fraglength)] 
    return slidingFragments

def editSeq(unsorted, dnaFragments):
    counts = sorted(unsorted)
    
    greatestOccurence = unsorted.index(counts[len(counts)-1])
    return [greatestOccurence, counts[len(counts)-1], dnaFragments[greatestOccurence]]


def generateAltSequence(codons,fragment, i):
    dnaSeq = ""
    for x in fragment: 
        if len(codons[x]) < i: 
            dnaSeq += codons[x][i]
        else:
            dnaSeq += codons[x][0]
    return dnaSeq 

def generateCodonUsage(codons, fragment): 
    generatedSequences = []
    for x in fragment:
        generatedSequences.append(generateAltSequence(codons, fragment, ))

    return generatedSequences

def adjSeq(fragmentedSeq, pattern):
    output = "" 
    for x in fragmentedSeq: 
        if x == pattern[2]:
            output += aatoDNA(x, 1)
        else: 
            output += aatoDNA(x, 0)

    return output

dnaSeq = aatoDNA(plist) 
print(dnaSeq)

dnaFragments = slidingFragmentation(p)
print(dnaFragments)

dnaFragmentCount = fragmentOccurence(dnaFragments, p)
print(dnaFragmentCount)


rep = editSeq(dnaFragmentCount, dnaFragments)

print(rep)
print(len(ecoli_codon_dict["A"]))    

fragmentationPattern = []
fragmentation = "" 

for x in p:
    fragmentation += x
    if rep[2] in fragmentation:
        fragmentationPattern.append(fragmentation[:-3])
        fragmentationPattern.append(fragmentation[-3:])
        fragmentation = "" 

print(fragmentationPattern)

print(adjSeq(fragmentationPattern, rep))