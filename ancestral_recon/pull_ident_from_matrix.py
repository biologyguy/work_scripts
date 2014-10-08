#!/usr/bin/python3
import sys


def pull_ident_from_matrix(matrix_file):
    top_line = matrix_file.readline()
    ids_list = top_line.strip().split(",")
    ids_list.pop(0)
    matrix_dict = {}
    
    for x in ids_list:
        matrix_dict[x] = {}
    
    for next_line in matrix_file:
        identities = next_line.strip().split(",")    
        seq_id = identities.pop(0)
        
        counter = 0
        for x in identities:
            matrix_dict[seq_id][ids_list[counter]] = x
            counter += 1
    
    return matrix_dict


matrix_file = open("seq_files/10000seq_subgroup_real_alignment.csv", "r")

matrix_dict = pull_ident_from_matrix(matrix_file)

matrix_file.close()

file_output = open("seq_files/real_alignment_IDs_for_histogram.txt", "w")
file_output.write("Comparison\t%ID\n")
complete_key = []

for key in matrix_dict:
    
    if key in complete_key:
        continue
    
    for d2key in matrix_dict[key]:
        if d2key == key:
            continue
        
        file_output.write(key + "->" + d2key + "\t" + matrix_dict[key][d2key] + "\n")
        
    complete_key.append(key)

file_output.close()
