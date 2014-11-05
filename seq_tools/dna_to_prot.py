#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Nov 5 2014 

"""
DESCRIPTION OF PROGRAM
"""

import argparse
import os
import shutil
import MyFuncs
from Bio import SeqIO
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(prog="dna_to_prot", description="Simple script that takes a DNA file and converts each to protein in reading frame #1")

parser.add_argument("in_file", help="Location of DNA fasta file", action="store")
parser.add_argument("-o", "--out_file", help="Send the protein seq to a new file. Otherwise it goes to std_out", action="store")
parser.add_argument("-f", "--format", help="If I get motivated, I'll let the in and out format be something other than fasta", type=str, choices=["fasta"], default="fasta")
parser.add_argument("-sd", "--strip_description", help="Clean up the crap stored in the description tag", action="store_true")
parser.add_argument("-ow", "--over_write", help="Replace a file that already exists with new out file", action="store_true")

in_args = parser.parse_args()

in_file = os.path.abspath(in_args.in_file)

prot_seqs = []
with open(in_file, "r") as ifile:
    dna_seqs = SeqIO.parse(ifile, "fasta")
    for seq in dna_seqs:
        if in_args.strip_description:
            seq.description = ""
        seq.alphabet = IUPAC.protein
        seq.seq = seq.seq.translate()
        prot_seqs.append(seq)

tmp_file = MyFuncs.TempFile()
with open(tmp_file.file, "w") as ofile:
    SeqIO.write(prot_seqs, ofile, "fasta")

if not in_args.out_file:
    with open(tmp_file.file, "r") as ifile:
        print(ifile.read())
else:
    out_file = os.path.abspath(in_args.out_file)
    if os.path.exists(out_file) and not in_args.over_write:
        print("Error: The outfile you've specified already exists. Use the -ow flag if you want to over-write it.")
    else:
        shutil.move(tmp_file.file, out_file)