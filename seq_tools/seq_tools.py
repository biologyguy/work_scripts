#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 

"""
DESCRIPTION OF PROGRAM
"""

import sys
import os
import re
from random import sample
from math import ceil
from Bio import SeqIO
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


# Apply DNA features to protein sequences
def map_features_dna2prot(dna_seqs, prot_seqs, quiet=False):  # Input as SeqIO.to_dict objects.
    new_seqs = {}
    for seq_id in dna_seqs:
        if seq_id not in prot_seqs:
            if not quiet:
                print("Warning: %s is in cDNA file, but not protein file" % seq_id)
            continue

        new_seqs[seq_id] = prot_seqs[seq_id]

        for feature in dna_seqs[seq_id].features:
            start = ceil(feature.location.start / 3)
            end = ceil(feature.location.end / 3)
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            prot_seqs[seq_id].features.append(new_feature)

    for seq_id in prot_seqs:
        if seq_id not in dna_seqs:
            if not quiet:
                print("Warning: %s is in protein file, but not the cDNA file" % seq_id)

    return new_seqs


# Apply DNA features to protein sequences
def map_features_prot2dna(prot_seqs, dna_seqs, quiet=False):  # Input as SeqIO.to_dict objects.
    new_seqs = {}
    for seq_id in prot_seqs:
        if seq_id not in dna_seqs:
            if not quiet:
                print("Warning: %s is in protein file, but not cDNA file" % seq_id)
            continue

        new_seqs[seq_id] = dna_seqs[seq_id]

        for feature in prot_seqs[seq_id].features:
            start = feature.location.start * 3 - 2
            end = feature.location.end * 3
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            dna_seqs[seq_id].features.append(new_feature)

    for seq_id in dna_seqs:
        if seq_id not in prot_seqs:
            if not quiet:
                print("Warning: %s is in cDNA file, but not protein file" % seq_id)

    return new_seqs


# Merge feature lists
def combine_features(seqs1, seqs2, quiet=False):  # Input as SeqIO.to_dict objects.
    # make sure that we're comparing apples to apples across all sequences (i.e., same alphabet)
    reference_alphabet = sample(seqs1.items(), 1)[0][1].seq.alphabet
    for seq_id in seqs1:
        if type(seqs1[seq_id].seq.alphabet) != type(reference_alphabet):
            print("You have mixed multiple alphabets into your sequences. Make sure everything is the same.")
            print("\t%s in first set" % seq_id)
            print("\tOffending alphabet: %s" % seqs1[seq_id].seq.alphabet)
            print("\tReference alphabet: %s" % reference_alphabet)
            sys.exit()

    for seq_id in seqs2:
        if type(seqs2[seq_id].seq.alphabet) != type(reference_alphabet):
            print("You have mixed multiple alphabets into your sequences. Make sure everything is the same.")
            print("\t%s in second set" % seq_id)
            print("\tOffending alphabet: %s" % seqs2[seq_id].seq.alphabet)
            print("\tReference alphabet: %s" % reference_alphabet)
            sys.exit()

    new_seqs = {}
    for seq_id in seqs1:
        if seq_id in seqs2:
            for feature in seqs2[seq_id].features:
                seqs1[seq_id].features.append(feature)
        else:
            if not quiet:
                print("Warning: %s is only in the first set of sequences" % seq_id)

        new_seqs[seq_id] = seqs1[seq_id]

    for seq_id in seqs2:
        if seq_id not in seqs1:
            if not quiet:
                print("Warning: %s is only in the first set of sequences" % seq_id)
            new_seqs[seq_id] = seqs2[seq_id]

    return new_seqs