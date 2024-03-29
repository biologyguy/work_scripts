#!/usr/bin/python3
from Bio import SeqIO
from math import sqrt 


def meanstdv(x): 
    n, mean, std = len(x), 0, 0 
    for a in x: 
        mean = mean + a 
    
    mean /= float(n) 
    
    for a in x: 
        std += (a - mean) ** 2 
    
    std = sqrt(std / float(n - 1)) 
    
    return mean, std


output = open("weird_sizes_removed.fasta", "w")
sequences = list(SeqIO.parse("dups_removed.fasta", "fasta"))

seq_length = []

for i in sequences:
    seq_length.append(len(i.seq))
    

AveStdev = meanstdv(seq_length)

seq_length = []
for i in sequences:
    lower = AveStdev[0] - AveStdev[1]
    upper = AveStdev[0] + AveStdev[1]
    length = len(i.seq)
    if lower < length < upper:
        SeqIO.write(i, output, 'fasta')
        seq_length.append(len(i.seq))
        
print(len(seq_length))
output.close()
