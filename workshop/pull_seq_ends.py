#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Oct 9 2014 

"""
Simple script to grab the front or end of every sequence in a fasta file and print a new fasta to stdout.
"""

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(prog='pull_seq_ends', description='')

parser.add_argument('seq_file', help='FASTA file', action='store')
parser.add_argument('amount', help='How much of the end of each sequence to you want', type=int, action='store')
parser.add_argument('which_end', help='beginning of sequences, or end?', choices=['front', 'end'])

in_args = parser.parse_args()
new_seqs = []
with open(in_args.seq_file, "r") as ifile:
    seqs = SeqIO.parse(ifile, "fasta")
    for seq in seqs:
        if in_args.which_end == 'front':
            seq.seq = seq.seq[:in_args.amount]
        else:
            seq.seq = seq.seq[-1 * in_args.amount:]
        print(seq)



