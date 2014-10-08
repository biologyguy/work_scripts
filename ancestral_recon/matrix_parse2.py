#!/usr/bin/python
from Bio import SeqIO 
import numpy as np
import sys
import argparse
import pdb
import os
from subprocess import Popen, PIPE, STDOUT

parser = argparse.ArgumentParser(prog="matrixParse.py", description="")
parser.add_argument('-g', '--gen_file', help='fasta file containing the seq_evolver generation of interest.')
parser.add_argument('-a', '--matrix_file', help='csv file with the matrix genious outputs from an alignment.')
parser.add_argument('-o', '--output_dir', help='Where to with the output?')
incoming_args = parser.parse_args()


def check_subgroup_lists(threshold_list, subgroup_lists):  # returns the index of a match in subgroup_lists, or false
    for thresh_id in threshold_list:
        subgroup_counter = 0
        for subgroup_list in subgroup_lists:
            if thresh_id in subgroup_list:
                return subgroup_counter
            subgroup_counter += 1
    return "append"


def merge_subgroup(subgroup_lists, idx, threshold_list):
    merged_subgroup = subgroup_lists[idx]
    for next_id in threshold_list:
        if next_id in merged_subgroup:
            continue
        else:
            merged_subgroup.append(next_id)
    
    subgroup_lists[idx] = merged_subgroup
    return subgroup_lists

matrix_file = open(incoming_args.matrix_file, "r")
top_line = matrix_file.readline()
ids_list = top_line.strip().split(",")
ids_list.pop(0)

matrix_dict = {}

print("Building matrix dictionary")
for x in ids_list:
    matrix_dict[x] = {}

percent_counter = 1.0
matrix_size = len(matrix_dict)

for next_line in matrix_file:
    identities = next_line.strip().split(",")    
    seq_id = identities.pop(0)
    
    percent_completion = (percent_counter / matrix_size) * 100
    sys.stdout.write("\rProcessing ID: " + seq_id + " (%.1f" % percent_completion + "%)")
    sys.stdout.flush()
    percent_counter += 1
    
    counter = 0
    for x in identities:
        matrix_dict[seq_id][ids_list[counter]] = x
        counter += 1
    

matrix_file.close()
print("\nIdentifying sequences with >35% identity")
#Go through all the rows and push the records with >35% Ident into a new dictionary
threshold_lists = {}
counter = 1.0
for y_axis in matrix_dict:
    percent_completion = counter / matrix_size * 100
    sys.stdout.write("\rProcessing ID: " + y_axis + " (%.1f" % percent_completion + "%)")
    sys.stdout.flush()
    if y_axis == "":
        print("Woops, where did that '' form from in y_axis?")
        sys.exit()
    
    threshold_lists[y_axis] = [y_axis]    
    for x_axis in matrix_dict[y_axis]:
        if x_axis == y_axis:
            continue
        
        if x_axis == "":
            print("Woops, where did that '' form from in x_axis?", y_axis)
            sys.exit()
        
        if float(matrix_dict[y_axis][x_axis]) > 35:
            threshold_lists[y_axis].append(x_axis)
    counter += 1


#Next I need to combine rows with overlap
print("\nCombining sequences into subgroups")
subgroup_lists = []
list_size = len(threshold_lists)
counter = 1.0
for row_id in threshold_lists:
    percent_completion = counter / list_size * 100
    sys.stdout.write("\rProcessing ID: " + row_id + " (%.1f" % percent_completion + "%)")
    sys.stdout.flush()
    counter += 1
    
    check_subgroups = check_subgroup_lists(threshold_lists[row_id], subgroup_lists)
    if check_subgroups == "append":
        subgroup_lists.append(threshold_lists[row_id])
    else:
        subgroup_lists = merge_subgroup(subgroup_lists, check_subgroups, threshold_lists[row_id])


#finally, get the sequences and write out the new subgroups to files. Also run PAGAN on all subgroups with 5+ sequences in it
print("\nAnalyzing " + str(len(subgroup_lists)) + " subgroup files")
seq_handle = SeqIO.parse(incoming_args.gen_file, "fasta")
seq_dict = SeqIO.to_dict(seq_handle)
sub_file_counter = 0
for subgroup_list in subgroup_lists:
    outfile_name = incoming_args.gen_file.split("/")[-1].split(".")[0] + "_sub" + str(sub_file_counter + 1) + ".fasta"
    sys.stdout.write(outfile_name)
    sys.stdout.flush()
    outfile_name = incoming_args.output_dir + outfile_name
    group_file = open(outfile_name, "w")

    for next_seq in subgroup_list:
        SeqIO.write(seq_dict[next_seq], group_file, "fasta")
    group_file.close()
    
    if len(subgroup_list) > 4:
        if not os.path.isdir(os.path.join(incoming_args.output_dir, "pagan")):
            sys.stdout.write(" ---> Making pagan directory")
            sys.stdout.flush()
            os.mkdir(os.path.join(incoming_args.output_dir, "pagan"))
        
        sys.stdout.write(" ---> executing PAGAN")
        sys.stdout.flush()
        
        event = Popen('pagan --seqfile ' + outfile_name + ' --outfile ' + incoming_args.output_dir + 'pagan/subgroup' + str(sub_file_counter + 1) + ' --output-ancestors', shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        event.communicate()
        print("")
    else:
        sys.stdout.write(" ---> only " + str(len(subgroup_list)) + " sequences")
        sys.stdout.flush()
    sub_file_counter += 1

print("\nDone")
