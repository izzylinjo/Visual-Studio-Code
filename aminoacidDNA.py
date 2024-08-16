
p = "MIMAAALLLAAALLLAAALLLAKLAKLAKLAALLAALLCL*"
#p2 = "MIMAAALAAALMMLAMMLAPPPLPPPLPPPPPPLLLLLL"


print(p.count("AAA"))


def aminoAcidFragmentation(aaseq, fraglength = 3, move = 0): 
    fragments = [(p[i: i+fraglength]) for i in range(move, len(p), fraglength)] 
    for i in range(0, len(fragments)): 
        if fragments[i] == fragments[i-1]:
            print("--------------------------")
            print("Replicate Detected")
            print(fragments[i])
            print("--------------------------")
        if len(fragments[i]) != len(fragments[i-1]):
            print("Fragment Not Equal Length")   

    return fragments

def shiftFragmentation(aaFragment, move=0):
    move  

def countOccurence(aaStr, fragment): 
    hello = 0 



for x in range(3,int(len(p)/2)+1):
   # print(aminoAcidFragmentation(p, fraglength=x))
    HELLO = 0 







