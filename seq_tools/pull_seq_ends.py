#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Oct 9 2014 

"""
Simple script to grab the front or end of every sequence in a fasta file and print a new fasta to stdout.
"""

import argparse
from Bio import SeqIO
import os
import MyFuncs
import shutil
import re

parser = argparse.ArgumentParser(prog='pull_seq_ends', description='')

parser.add_argument('seq_file', help='FASTA file', action='store')
parser.add_argument('amount', help='How much of the end of each sequence to you want', type=int, action='store')
parser.add_argument('which_end', help='beginning of sequences, or end?', choices=['front', 'end'])
parser.add_argument('in_format', help='What is the format of your sequence?', choices=['fasta', 'phylip', 'phylip-relaxed'], default='fasta')
parser.add_argument('out_format', help='What format do you want output?', choices=['fasta', 'phylip', 'phylip-relaxed'], default='fasta')
parser.add_argument('-o', '--outfile', action='store')

in_args = parser.parse_args()
new_seqs = []
with open(in_args.seq_file, "r") as ifile:
    seqs = SeqIO.parse(ifile, in_args.in_format)
    for seq in seqs:
        if in_args.which_end == 'front':
            seq.seq = seq.seq[:in_args.amount]

        else:
            seq.seq = seq.seq[-1 * in_args.amount:]

        new_seqs.append(seq)

temp_file = MyFuncs.TempFile()
with open(temp_file.file, "w") as ofile:
    SeqIO.write(new_seqs, ofile, in_args.out_format)

# Need to fix the phylip file...
if in_args.out_format in ["phylip", "phylip-relaxed"]:
    with open(temp_file.file, "r") as ifile:
        num_seqs, length = ifile.readline().strip().split(" ")
        chunks = ifile.read().split("\n\n")
        joined_seqs = []
        for i in range(int(num_seqs)):
            joined_seqs.append("")

        for i in range(len(chunks)):
            seqs = chunks[i].strip().split("\n")
            for j in range(len(seqs)):
                seq = seqs[j].strip()
                seq = re.sub(r"([ATCG-]) ([ATCG-])", r"\1\2", seq)
                joined_seqs[j] += seq

        new_phylip = " %s %s\n" % (num_seqs, length)
        for seq in joined_seqs:
            new_phylip += "%s\n" % seq

    with open(temp_file.file, "w") as ofile:
        ofile.write(new_phylip)

if in_args.outfile:
    shutil.copyfile(temp_file.file, os.path.abspath(in_args.outfile))
    print("New file written: %s" % os.path.abspath(in_args.outfile))

else:
    with open(temp_file.file, "r") as ifile:
        print(ifile.read())