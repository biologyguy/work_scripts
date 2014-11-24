#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

import argparse
import sys
import re
from Bio import AlignIO
from Bio.SubsMat import MatrixInfo


parser = argparse.ArgumentParser(prog="msa_parse", description="Takes a multiple sequence alignment, and "
                                                               "outputs the %ID and %SIM for all vs. all "
                                                               "(or subset if -s flag is used)")
parser.add_argument('in_file', help='Must be in fasta format', action='store')
parser.add_argument('-s', '--search_string',
                    help='groups a subset of sequences together, and compares them against everyone else',
                    action="store", default=False)
parser.add_argument('-i', '--internal', help='Only compares the sequences in -s', action="store_true", default=False)
parser.add_argument('-n', '--not_internal', help='Only compares the sequences not in -s',
                    action="store_true", default=False)

in_args = parser.parse_args()


def convert_sub_matrix(sub_matrix):
        """
        this is taking a MatrixInfo object, and converting it into a 2D dictionary instead of a 1D dictionary, stripping
        out residues B and Z, and adding in gap penalty. It will break if some other type of data is fed into it.
        """
        converted_sub_matrix = {"A": {}, "C": {}, "D": {}, "E": {}, "F": {}, "G": {}, "H": {}, "I": {}, "K": {},
                                "L": {}, "M": {}, "N": {}, "P": {}, "Q": {}, "R": {}, "S": {}, "T": {}, "V": {},
                                "W": {}, "Y": {}, "X": {}, "-": {}}
        for next_aa in converted_sub_matrix:
            if next_aa == "-":
                converted_sub_matrix[next_aa] = {"A": (-4), "C": (-4), "D": (-4), "E": (-4), "F": (-4), "G": (-4),
                                                 "H": (-4), "I": (-4), "K": (-4), "L": (-4), "M": (-4), "N": (-4),
                                                 "P": (-4), "Q": (-4), "R": (-4), "S": (-4), "T": (-4), "V": (-4),
                                                 "W": (-4), "Y": (-4), "X": (-4), "-": (-4)}
            elif next_aa == "X":
                converted_sub_matrix[next_aa] = {"A": (-1), "C": (-1), "D": (-1), "E": (-1), "F": (-1), "G": (-1),
                                                 "H": (-1), "I": (-1), "K": (-1), "L": (-1), "M": (-1), "N": (-1),
                                                 "P": (-1), "Q": (-1), "R": (-1), "S": (-1), "T": (-1), "V": (-1),
                                                 "W": (-1), "Y": (-1), "X": (-1), "-": (-4)}
            else:
                converted_sub_matrix[next_aa] = {"A": 0, "C": 0, "D": 0, "E": 0, "F": 0, "G": 0, "H": 0, "I": 0,
                                                 "K": 0, "L": 0, "M": 0, "N": 0, "P": 0, "Q": 0, "R": 0, "S": 0,
                                                 "T": 0, "V": 0, "W": 0, "Y": 0, "X": (-1), "-": (-4)}
                   
        for k in sub_matrix:
            x = k[0]
            y = k[1]
            
            if x in ['B', 'Z', 'X'] or y in ['B', 'Z', 'X']:
                continue
            
            converted_sub_matrix[x][y] = sub_matrix[k]
            converted_sub_matrix[y][x] = sub_matrix[k]
               
        return converted_sub_matrix
 
 
def comp(aln1, aln2, in_matrix):
    seq1 = list(aln1.seq)
    seq2 = list(aln2.seq)
    
    output = {'id': 0.0, "sim": 0.0, "len": 0.0}
    for index in range(len(seq1)):
        if seq1[index] == '-' and seq2[index] == '-':
            continue
        elif seq1[index] == '-' or seq2[index] == '-':
            output['len'] += 1
            continue
        elif seq1[index] == seq2[index]:
            output['id'] += 1.0
            output['sim'] += 1.0
            output['len'] += 1
            continue
        elif int(in_matrix[seq1[index]][seq2[index]]) >= 0:
            output['sim'] += 1.0
            output['len'] += 1
            continue
        else:
            output['len'] += 1.0
        
    output['sim'] /= output['len']
    output['id'] /= output['len']
    
    return output

with open(in_args.in_file, 'r') as alignment_file:
    align_obj = AlignIO.read(alignment_file, 'fasta')

matrix = convert_sub_matrix(MatrixInfo.blosum62)
output_dir = {}  # [{(Name1, Name2): {'id':%ID, 'sim':%SIM, 'len':LEN}]

group1 = []
clip_group1 = False
group2 = []

if (in_args.internal or in_args.not_internal) and not in_args.search_string:
    sys.exit("You can't set the -i or -n flag without setting -s as well...")

if in_args.internal and in_args.not_internal:
    sys.exit("You can't set both -i and -n flags at the same time...")
    
if in_args.internal:  
    for i in align_obj:
        if re.search(in_args.search_string, i.id):
            group1.append(i)
            group2.append(i)
            clip_group1 = True

elif in_args.not_internal:  
    for i in align_obj:
        if not re.search(in_args.search_string, i.id):
            group1.append(i)
            group2.append(i)
            clip_group1 = True
                
elif in_args.search_string:
    for i in align_obj:
        if re.search(in_args.search_string, i.id):
            group1.append(i)
        else:
            group2.append(i)

else:
    group1 = align_obj
    group2 = align_obj
    clip_group1 = True
    
number_comparisons = len(group1)
print("Rec1\tRec2\t%ID\t%Sim\tLength")

for i in group1:
    if clip_group1:
        del group2[0]
    for j in group2:
        if i.id == j.id:
            continue
        comparison = comp(i, j, matrix)
        output_dir[(i.id, j.id)] = comparison
        
        print(i.id, "\t", j.id, "\t", round(comparison['id'], 3), "\t", round(comparison['sim'], 3),
              "\t", comparison['len'])