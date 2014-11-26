#!/usr/bin/python3
from Bio import SeqIO

sequences = SeqIO.parse("genome_protein_sequences/Acropora_digitifera_v1.0.1.prot.fa", "fasta")
output = open("Acropora_digitifera_v1.0.1.prot.fa", "w")

ids_dict = {}

for i in sequences:
    if i.id in ids_dict:
        ids_dict[i.id] += 1
    else:
        ids_dict[i.id] = 1
    
    i.id += "_" + str(ids_dict[i.id])
    SeqIO.write(i, output,'fasta')