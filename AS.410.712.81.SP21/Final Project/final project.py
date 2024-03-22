#!/usr/bin/env python3


import re

def reverse_complement(oligo_seq):
    oligo_join=(','.join(oligo_seq))
    oligo_seq_split=oligo_join.split(',')
    oligo_seq_split.reverse()
    n=len(oligo_seq)

    for x in range(0, n):
        if oligo_seq_split[x]=="C":
            oligo_seq_split[x]="G"
        elif oligo_seq_split[x]=="G":
            oligo_seq_split[x]="C"
        elif oligo_seq_split[x]=="A":
            oligo_seq_split[x]="T"
        elif oligo_seq_split[x]=="T":
            oligo_seq_split[x]="A"
    str = ""
    oligo_seq_rev = str.join( oligo_seq_split )
    return oligo_seq_rev

def GC(oligo_seq):
    G = oligo_seq.count("G")
    C = oligo_seq.count("C")
    GC_content=((G+C)/len(oligo_seq))*100
    return GC_content

def Tm(oligo_seq):
    G = oligo_seq.count("G")
    C = oligo_seq.count("C")
    A = oligo_seq.count("A")
    T = oligo_seq.count("T")
    if len(oligo_seq)> 14:
        Tm= 64.9 +41*(G+C-16.4)/(A+T+G+C)
    else:
        Tm= (A+T)*2 + (G+C)*4
    return Tm



oligo_seq = input("Enter an oligonucleotide sequence: ") 
if re.search(r"[^ATGC]", oligo_seq):
    print("An invalid character was found!")

else:
    print("Sequence (5' -> 3'): " + oligo_seq)
    print("Reverse complement (5' -> 3'): " +reverse_complement(oligo_seq))
    print("Length: " +str(len(oligo_seq)))
    print("GC content[%]: " + str(GC(oligo_seq))+"%")
    print("Melting temperature[°C]: " +str('%.1f' % Tm(oligo_seq)) + "°C")




