#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 

"""
DESCRIPTION OF PROGRAM
"""

import sys, os, re, shutil, MyFuncs
from math import ceil
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation


def guess_alphabet(sequence):  # Can be fasta file or raw, does not handle ambigious dna
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
    sequence = re.sub("[0-9\s\-\*]", "", sequence)
    sequence = sequence.upper()
    return sequence


def map_features_dna2prot(dna_file, prot_file, outfile=None):  # dna_file in gb, prot_file in fasta, outfile in gb
    with open(dna_file, "r") as ifile:
        cdna_seqs = SeqIO.to_dict(SeqIO.parse(ifile, "gb", IUPAC.ambiguous_dna))

    with open(prot_file, "r") as ifile:
        prot_seqs = SeqIO.to_dict(SeqIO.parse(ifile, "fasta", IUPAC.protein))

    new_seqs = {}
    for seq_id in cdna_seqs:
        if seq_id not in prot_seqs:
            print("Warning: %s is in cDNA file, but not protein file" % seq_id)
            continue

        new_seqs[seq_id] = prot_seqs[seq_id]

        for feature in cdna_seqs[seq_id].features:
            start = ceil(feature.location.start / 3)
            end = ceil(feature.location.end / 3)
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            prot_seqs[seq_id].features.append(new_feature)

    for seq_id in prot_seqs:
        if seq_id not in cdna_seqs:
            print("Warning: %s is in protein file, but not the cDNA file" % seq_id)

    new_seqs = [k for k in new_seqs.values()]

    if outfile:
        with open(outfile, "w") as ofile:
            SeqIO.write(new_seqs, ofile, "gb")
    else:
        for seq in new_seqs:
            print(seq.format("gb"))


def map_features_prot2dna(prot_file, dna_file, outfile=None):  # prot_file in gb, dna_file in fasta, outfile in gb
    with open(prot_file, "r") as ifile:
        prot_seqs = SeqIO.to_dict(SeqIO.parse(ifile, "gb", IUPAC.protein))

    with open(dna_file, "r") as ifile:
        cdna_seqs = SeqIO.to_dict(SeqIO.parse(ifile, "fasta", IUPAC.unambiguous_dna))

    new_seqs = {}
    for seq_id in prot_seqs:
        if seq_id not in cdna_seqs:
            print("Warning: %s is in protein file, but not cDNA file" % seq_id)
            continue

        new_seqs[seq_id] = cdna_seqs[seq_id]

        for feature in prot_seqs[seq_id].features:
            start = feature.location.start * 3 - 2
            end = feature.location.end * 3
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            cdna_seqs[seq_id].features.append(new_feature)

    for seq_id in cdna_seqs:
        if seq_id not in prot_seqs:
            print("Warning: %s is in cDNA file, but not protein file" % seq_id)

    new_seqs = [k for k in new_seqs.values()]

    if outfile:
        with open(outfile, "w") as ofile:
            SeqIO.write(new_seqs, ofile, "gb")
    else:
        for seq in new_seqs:
            print(seq.format("gb"))