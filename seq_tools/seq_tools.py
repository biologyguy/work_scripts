#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 

"""
Collection of functions that do funs stuff with sequences. Pull them into a script, or run as a commandline tool.
"""
import pdb
import sys
import os
import re
from random import sample
from math import ceil
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import AlignIO
from MyFuncs import TempFile


def guess_alphabet(sequence):  # Can be fasta file or raw, does not handle ambigious dna
    if os.path.isfile(sequence):
        seq_format = guess_format(sequence)
        if not seq_format:
            sys.exit("Error: could not determine the format of your input sequence file.")
        with open(sequence, "r") as infile:
            sequences = SeqIO.parse(infile, seq_format)
            sequence = ""
            for next_seq in sequences:
                if len(sequence) > 1000:
                    break
                sequence += str(next_seq.seq)

    sequence = clean_seq(sequence)
    sequence = re.sub("[NX]", "", sequence)

    percent_dna = float(sequence.count("A") + sequence.count("G") +
                        sequence.count("T") + sequence.count("C")) / float(len(sequence))
    if percent_dna > 0.95:
        return "nucl"
    else:
        return "prot"


def clean_seq(sequence):  # fasta file or raw
    """remove fasta headers, numbers, and whitespace from sequence string"""
    if os.path.isfile(sequence):
        with open(sequence, "r") as infile:
            sequence = SeqIO.read(infile, "fasta")
            sequence = str(sequence.seq)

    sequence = re.sub(">.*", "", sequence)
    sequence = re.sub("[0-9\s\-\*]", "", sequence)
    sequence = sequence.upper()
    return sequence


def concat_seqs(sequences):
    print("not implemented")


def combine_files(file_paths, mix=False):  # Combine the sequences from a bunch of files into a single list
    new_seq_list = []
    dna_or_prot = None
    for next_file in file_paths:
        in_format = guess_format(next_file)
        with open(os.path.abspath(next_file), "r") as infile:
            sequences = list(SeqIO.parse(infile, in_format))
            if not dna_or_prot:
                dna_or_prot = guess_alphabet(str(sequences[0].seq))
            for next_seq in sequences:
                alpha = guess_alphabet(str(next_seq.seq))
                if dna_or_prot != alpha and not mix:
                    sys.exit("Error: It looks like you are trying to mix DNA and Protein files. If you really do want"
                             "to do that, set mix=True")

                if alpha == "nucl":
                    next_seq.seq.alphabet = IUPAC.ambiguous_dna

                if alpha == "prot":
                    next_seq.seq.alphabet = IUPAC.protein
                new_seq_list.append(next_seq)

    return new_seq_list


def guess_format(in_file):
    # currently just looks at file extension, but might want to get more fancy
    if not os.path.isfile(in_file):
        return None

    in_file = in_file.split(".")
    if in_file[-1] in ["fa", "fas", "fasta"]:
        return "fasta"
    if in_file[-1] in ["gb", "genbank"]:
        return "gb"
    if in_file[-1] in ["nex", "nxs", "nexus"]:
        return "nexus"
    if in_file[-1] in ["phy", "ph", "phylip"]:
        return "phylip-relaxed"
    if in_file[-1] in ["sto", "pfam", "stockholm"]:
        return "stockholm"
    return None


# Apply DNA features to protein sequences
def map_features_dna2prot(dna_seqs, prot_seqs, quiet=False):  # Input as SeqIO.to_dict objects.
    _new_seqs = {}
    for _seq_id in dna_seqs:
        if _seq_id not in prot_seqs:
            if not quiet:
                print("Warning: %s is in cDNA file, but not protein file" % _seq_id, file=sys.stderr)
            continue

        _new_seqs[_seq_id] = prot_seqs[_seq_id]

        for feature in dna_seqs[_seq_id].features:
            start = ceil(feature.location.start / 3)
            end = ceil(feature.location.end / 3)
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            prot_seqs[_seq_id].features.append(new_feature)

    for _seq_id in prot_seqs:
        if _seq_id not in dna_seqs:
            if not quiet:
                print("Warning: %s is in protein file, but not the cDNA file" % _seq_id, file=sys.stderr)

    return _new_seqs


# Apply DNA features to protein sequences
def map_features_prot2dna(prot_seqs, dna_seqs, quiet=False):  # Input as SeqIO.to_dict objects.
    _new_seqs = {}
    for _seq_id in prot_seqs:
        if _seq_id not in dna_seqs:
            if not quiet:
                print("Warning: %s is in protein file, but not cDNA file" % _seq_id)
            continue

        _new_seqs[_seq_id] = dna_seqs[_seq_id]

        for feature in prot_seqs[_seq_id].features:
            start = feature.location.start * 3 - 2
            end = feature.location.end * 3
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            dna_seqs[_seq_id].features.append(new_feature)

    for _seq_id in dna_seqs:
        if _seq_id not in prot_seqs:
            if not quiet:
                print("Warning: %s is in cDNA file, but not protein file" % _seq_id)

    return _new_seqs


# Merge feature lists
def combine_features(seqs1, seqs2, quiet=False):  # Input as SeqIO.to_dict objects.
    # make sure that we're comparing apples to apples across all sequences (i.e., same alphabet)
    reference_alphabet = sample(seqs1.items(), 1)[0][1].seq.alphabet
    for _seq_id in seqs1:
        if type(seqs1[_seq_id].seq.alphabet) != type(reference_alphabet):
            print("You have mixed multiple alphabets into your sequences. Make sure everything is the same.")
            print("\t%s in first set" % _seq_id)
            print("\tOffending alphabet: %s" % seqs1[_seq_id].seq.alphabet)
            print("\tReference alphabet: %s" % reference_alphabet)
            sys.exit()

    for _seq_id in seqs2:
        if type(seqs2[_seq_id].seq.alphabet) != type(reference_alphabet):
            print("You have mixed multiple alphabets into your sequences. Make sure everything is the same.")
            print("\t%s in second set" % _seq_id)
            print("\tOffending alphabet: %s" % seqs2[_seq_id].seq.alphabet)
            print("\tReference alphabet: %s" % reference_alphabet)
            sys.exit()

    _new_seqs = {}
    for _seq_id in seqs1:
        if _seq_id in seqs2:
            for feature in seqs2[_seq_id].features:
                seqs1[_seq_id].features.append(feature)
        else:
            if not quiet:
                print("Warning: %s is only in the first set of sequences" % _seq_id)

        _new_seqs[_seq_id] = seqs1[_seq_id]

    for _seq_id in seqs2:
        if _seq_id not in seqs1:
            if not quiet:
                print("Warning: %s is only in the first set of sequences" % _seq_id)
            _new_seqs[_seq_id] = seqs2[_seq_id]

    return _new_seqs


def print_format(sequences, file_format):

    for sequence in sequences:
        print(sequence.format(file_format))


if __name__ == '__main__':
    import argparse
    from Bio.Alphabet import IUPAC

    parser = argparse.ArgumentParser(prog="seq_tools.py", description="Commandline wrapper for all the fun functions in"
                                                                      "this file. Play with your sequences!")

    parser.add_argument('-ga', '--guess_alphabet', action='store')
    parser.add_argument('-gf', '--guess_format', action='store')
    parser.add_argument('-cs', '--clean_seq', action='store')
    parser.add_argument('-fd2p', '--map_features_dna2prot', action='store', nargs=2)
    parser.add_argument('-fp2d', '--map_features_prot2dna', action='store', nargs=2)
    parser.add_argument('-cf', '--combine_features', action='store', nargs=2)
    parser.add_argument('-cl', '--combine_files', action='store', nargs="+", help="Arguments need to be supplied "
                                                                                  "<format> <files ... > <True|False "
                                                                                  "(for mix, optional)>")
    parser.add_argument('-sf', '--screw_formats', action='store', nargs=2, help="Arguments: <in_file> <out_format>")
    parser.add_argument('-sfa', '--screw_formats_align', action='store', nargs=2, help="Arguments: <in_file> <out_format>")

    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')

    in_args = parser.parse_args()

    # Screw formats
    if in_args.screw_formats:
        with open(os.path.abspath(in_args.screw_formats[0]), "r") as ifile:
            seqs = list(SeqIO.parse(ifile, guess_format(in_args.screw_formats[0])))

        for seq in seqs:
            print(seq.format(in_args.screw_formats[1]))

    # Screw formats align
    if in_args.screw_formats_align:
        tmp_file = TempFile()
        with open(os.path.abspath(in_args.screw_formats_align[0]), "r") as ifile:

            alignments = AlignIO.parse(ifile, guess_format(in_args.screw_formats_align[0]))
            with open(tmp_file.file, "w") as ofile:
                AlignIO.write(alignments, ofile, in_args.screw_formats_align[1])

        with open(tmp_file.file, "r") as ifile:
            print(ifile.read())

    # Guess alphabet
    if in_args.guess_alphabet:
        print(guess_alphabet(in_args.guess_alphabet))

    # Clean Seq
    if in_args.clean_seq:
        print(clean_seq(in_args.clean_seq))

    # Guess format
    if in_args.guess_format:
        print(guess_format(in_args.guess_format))

    # Map features from cDNA over to protein
    if in_args.map_features_dna2prot:
        dna, prot = in_args.map_features_dna2prot
        dna = os.path.abspath(dna)
        prot = os.path.abspath(prot)

        if guess_format(dna) != "gb":
            sys.exit("Error: first argument must be genbank format")
        if guess_alphabet(dna) != "nucl":
            sys.exit("Error: first argument must be the dna sequence")
        if not guess_format(prot):
            sys.exit("Error: couldn't determine the format of your protein sequence")
        if guess_alphabet(prot) != "prot":
            sys.exit("Error: second argument must be the protein sequence")

        with open(dna, "r") as ifile:
            dna = SeqIO.to_dict(SeqIO.parse(ifile, "gb"))

        with open(prot, "r") as ifile:
            prot = SeqIO.to_dict(SeqIO.parse(ifile, guess_format(prot)))

        new_seqs = map_features_dna2prot(dna, prot)
        for seq_id in new_seqs:
            new_seqs[seq_id].seq.alphabet = IUPAC.protein
            print(new_seqs[seq_id].format("gb"))

    # Map features from protein over to cDNA
    if in_args.map_features_prot2dna:
        prot, dna = in_args.map_features_prot2dna
        dna = os.path.abspath(dna)
        prot = os.path.abspath(prot)

        if guess_format(prot) != "gb":
            sys.exit("Error: first argument must be genbank format")
        if guess_alphabet(prot) != "prot":
            sys.exit("Error: first argument must be the protein sequence")
        if not guess_format(dna):
            sys.exit("Error: couldn't determine the format of your dna sequence")
        if guess_alphabet(dna) != "nucl":
            sys.exit("Error: second argument must be the dna sequence")

        with open(prot, "r") as ifile:
            prot = SeqIO.to_dict(SeqIO.parse(ifile, "gb"))

        with open(dna, "r") as ifile:
            dna = SeqIO.to_dict(SeqIO.parse(ifile, guess_format(dna)))

        new_seqs = map_features_prot2dna(prot, dna)
        for seq_id in new_seqs:
            new_seqs[seq_id].seq.alphabet = IUPAC.ambiguous_dna
            print(new_seqs[seq_id].format("gb"))

    # Combine feature sets from two files into one
    if in_args.combine_features:
        file1, file2 = in_args.combine_features
        file1 = os.path.abspath(file1)
        file2 = os.path.abspath(file2)

        if guess_format(file1) != "gb" or guess_format(file2) != "gb":
            sys.exit("Error: please provide sequence files with the .gb extension")

        with open(file1, "r") as ifile:
            file1 = SeqIO.to_dict(SeqIO.parse(ifile, "gb"))

        with open(file2, "r") as ifile:
            file2 = SeqIO.to_dict(SeqIO.parse(ifile, "gb"))

        new_seqs = combine_features(file1, file2)
        for seq_id in new_seqs:
            new_seqs[seq_id].seq.alphabet = IUPAC.protein
            print(new_seqs[seq_id].format("gb"))

    # Combine group of files into one
    if in_args.combine_files:
        if in_args.combine_files[-1].upper() == "TRUE":
            mix_alphabet = True
            files = in_args.combine_files[1:-1]
        elif in_args.combine_files[-1].upper() == "FALSE":
            mix_alphabet = False
            files = in_args.combine_files[1:-1]
        else:
            mix_alphabet = False
            files = in_args.combine_files[1:]

        new_seqs = combine_files(files, mix_alphabet)
        for seq in new_seqs:
            print(seq.format(in_args.combine_files[0]))