#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 

"""
DESCRIPTION OF PROGRAM
"""

import sys, os, re, shutil, MyFuncs
from Bio import SeqIO


def guess_alphabet(sequence):  # Can be fasta file or raw
    if os.path.isfile(sequence):
        with open(sequence, "r") as ifile:
            sequences = SeqIO.parse(ifile, "fasta")
            sequence = ""
            for seq in sequences:
                if len(sequence) > 1000:
                    break
                sequence += seq.seq

    sequence = clean_seq(sequence)
    sequence = re.sub("[NX]", "", sequence)

    percent_dna = float(sequence.count("A") + sequence.count("G") + sequence.count("T") + sequence.count("C")) / float(len(sequence))
    if percent_dna > 0.95:
        return "nucl"
    else:
        return "prot"


def clean_seq(sequence):
    """remove fasta headers, numbers, and whitespace from sequence string"""
    sequence = re.sub(">.*", "", sequence)
    sequence = re.sub("[0-9\s]", "", sequence)
    sequence = sequence.upper()
    return sequence
