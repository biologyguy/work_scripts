#!/usr/bin/env python3

# Little script Bernie wrote to pull a singe fasta record out of a file given part or all of the heading info
from Bio import SeqIO
from sys import argv, exit
from os import path
if len(argv) != 3 or argv[1] == '-h' or argv[1] == '-help':
    print("Christy_script.py <fasta file> <substring>")
    print("Enter substring to search for in fasta headers.\nReturns all fasta sequences with substring in header")
    exit(0)

if path.isfile(argv[1]):
    fasta_file = argv[1]
    substring = argv[2]

else:
    fasta_file = argv[2]
    substring = argv[1]

for seq in SeqIO.parse(fasta_file, 'fasta'):
    if seq.description.find(substring) != -1:
        print(seq.format('fasta'))
exit(0)