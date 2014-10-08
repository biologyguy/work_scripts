#!/usr/bin/python
from Bio import SeqIO 
import argparse



matrix_file = open("seq_files/gen100_align_matrix.txt", "r")
top_line = matrix_file.readline()
ids_list = top_line.strip().split(",")
ids_list.pop(0)

seq_handle = SeqIO.parse("seq_files/gen100.fasta", "fasta")
seq_dict = SeqIO.to_dict(seq_handle)

sub_file_counter = 0
progression_count = 0

for next_line in matrix_file:
    identities = next_line.strip().split(",")
    
    if progression_count > 0:
        progression_count -= 1
        continue
    
    group_file = open("seq_files/gen100_sub" + str(sub_file_counter) + ".fasta", "w")
    sub_file_counter += 1
    
    seq_id = identities.pop(0)
    SeqIO.write(seq_dict[seq_id], group_file, "fasta")
    progression_count = 0
    index_counter = -1
    for next_value in identities:
        index_counter += 1
        if next_value == "":
            continue
        
        if float(next_value) > 35.0:
            progression_count += 1
            SeqIO.write(seq_dict[str(ids_list[index_counter])], group_file, "fasta")
        
        else:
            continue


matrix_file.close()