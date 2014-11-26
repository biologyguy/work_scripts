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
import string
from random import sample, choice
from math import ceil
from Bio import SeqIO, AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError


# ################################################ INTERNAL FUNCTIONS ################################################ #
def _set_alphabet(_sequence):  # update sequence alphabet in place
    alpha = guess_alphabet(str(_sequence))
    if alpha == "nucl":
        _sequence.alphabet = IUPAC.ambiguous_dna
    elif alpha == "prot":
        _sequence.alphabet = IUPAC.protein
    else:
        sys.exit("Error: Can't deterimine alphabet in _set_alphabet")
    return _sequence


def _sequence_list(sequence):  # Open a file and parse, or convert raw into a Seq object
    if os.path.isfile(sequence):
        _seq_format = guess_format(sequence)
        if not _seq_format:
            sys.exit("Error: could not determine the format of your input sequence file.")
        with open(sequence, "r") as infile:
            _sequences = list(SeqIO.parse(infile, _seq_format))
    else:
        # dna_or_prot = IUPAC.protein if guess_alphabet(sequence) == "prot" else IUPAC.ambiguous_dna
        _sequences = [Seq(sequence)]

    return _sequences
# #################################################################################################################### #


def guess_alphabet(sequence):  # Can be sequence file or raw, does not handle ambigious dna
    _sequences = _sequence_list(sequence)
    sequence = ""
    for next_seq in _sequences:
        if len(sequence) > 1000:
            break
        sequence += str(next_seq.seq)

    sequence = re.sub("[NX]", "", sequence)

    percent_dna = float(sequence.count("A") + sequence.count("G") +
                        sequence.count("T") + sequence.count("C")) / float(len(sequence))
    if percent_dna > 0.95:
        return "nucl"
    else:
        return "prot"


def translate_cds(_sequences):
    _output = []
    _sequences = _sequence_list(_sequences)

    for _seq in _sequences:
        try:
            _seq.seq = _seq.seq.translate(cds=True, to_stop=True)
        except TranslationError as e1:
            _seq.seq = Seq(str(_seq.seq)[:(len(str(_seq.seq)) - len(str(_seq.seq)) % 3)])
            try:
                _seq.seq = _seq.seq.translate()
                print("Warning: %s is not a standard CDS\t-->\t%s" % (_seq.id, e1), file=sys.stderr)
            except TranslationError as e2:
                print("Error: %s failed to translate\t-->\t%s" % (_seq.id, e2), file=sys.stderr)

        _seq.seq.alphabet = IUPAC.protein
        _output.append(_seq)
    return _output


def concat_seqs(_sequences):
    _sequences = clean_seq(_sequences)
    return "".join(_sequences)


def count_sequences(_sequences):
    _sequences = _sequence_list(_sequences)
    return len(_sequences)


def clean_seq(sequence):  # from file or raw
    """remove fasta headers, numbers, and whitespace from sequence strings"""
    _sequences = _sequence_list(sequence)
    _output = []
    for _seq in _sequences:
        _seq = re.sub(">.*", "", str(_seq.seq))
        _seq = re.sub("[0-9\s\-\*]", "", _seq)
        _seq = _seq.upper()
        _output.append(_seq)

    return _output  # returns a list of sequences in string format


def combine_files(file_paths, mix=False):  # Combine the sequences from a bunch of files into a single list
    new_seq_list = []
    dna_or_prot = None
    for next_file in file_paths:
        in_format = guess_format(next_file)
        with open(os.path.abspath(next_file), "r") as infile:
            _sequences = list(SeqIO.parse(infile, in_format))
            if not dna_or_prot:
                dna_or_prot = guess_alphabet(str(_sequences[0].seq))
            for next_seq in _sequences:
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
def map_features_prot2dna(prot_seqs, dna_seqs):  # Input as SeqIO.to_dict objects.
    _new_seqs = {}
    for _seq_id in prot_seqs:
        if _seq_id not in dna_seqs:
            print("Warning: %s is in protein file, but not cDNA file" % _seq_id, file=sys.stderr)
            continue

        _new_seqs[_seq_id] = dna_seqs[_seq_id]

        for feature in prot_seqs[_seq_id].features:
            start = feature.location.start * 3 - 2
            end = feature.location.end * 3
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            dna_seqs[_seq_id].features.append(new_feature)

    for _seq_id in dna_seqs:
        if _seq_id not in prot_seqs:
            print("Warning: %s is in cDNA file, but not protein file" % _seq_id, file=sys.stderr)

    return _new_seqs


# Merge feature lists
def combine_features(seqs1, seqs2):  # Input as SeqIO.to_dict objects.
    # make sure that we're comparing apples to apples across all sequences (i.e., same alphabet)
    reference_alphabet = sample(seqs1.items(), 1)[0][1].seq.alphabet
    for _seq_id in seqs1:
        if type(seqs1[_seq_id].seq.alphabet) != type(reference_alphabet):
            error_mes = "You have mixed multiple alphabets into your sequences. Make sure everything is the same.\n" \
                        "\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s" \
                        % (_seq_id, seqs1[_seq_id].seq.alphabet, reference_alphabet)
            print(error_mes, file=sys.stderr)
            return False

    for _seq_id in seqs2:
        if type(seqs2[_seq_id].seq.alphabet) != type(reference_alphabet):
            error_mes = "You have mixed multiple alphabets into your sequences. Make sure everything is the same.\n" \
                        "\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s" \
                        % (_seq_id, seqs2[_seq_id].seq.alphabet, reference_alphabet)
            print(error_mes, file=sys.stderr)
            return False

    _new_seqs = {}
    for _seq_id in seqs1:
        if _seq_id in seqs2:
            for feature in seqs2[_seq_id].features:
                seqs1[_seq_id].features.append(feature)
        else:
            print("Warning: %s is only in the first set of sequences" % _seq_id, file=sys.stderr)

        _new_seqs[_seq_id] = seqs1[_seq_id]

    for _seq_id in seqs2:
        if _seq_id not in seqs1:
            print("Warning: %s is only in the first set of sequences" % _seq_id, file=sys.stderr)
            _new_seqs[_seq_id] = seqs2[_seq_id]

    return _new_seqs


def screw_formats(_sequences, file_format):
    return _sequences.format(file_format)


def screw_formats_align(_alignments, _out_format):
    _output = ""
    if _out_format == "phylipi":
        if len(_alignments) > 1:
            print("Warning: the input file contained more than one alignment, but phylip can only handle one. "
                  "The topmost alignment is shown here.", file=sys.stderr)
        _seqs = list(_alignments[0])
        _output += " %s %s\n" % (len(_seqs), len(_seqs[0].seq))
        max_id_length = 0
        for _seq in _seqs:
            max_id_length = len(_seq.id) if len(_seq.id) > max_id_length else max_id_length

        for _seq in _seqs:
            _seq_id = _seq.id.ljust(max_id_length)
            _output += "%s %s\n" % (_seq_id, _seq.seq)
    else:
        for alignment in _alignments:
            _output += alignment.format(_out_format)

    return _output


def hash_seqeunce_ids(_sequences):
    hash_list = []
    seq_ids = []
    for i in range(len(_sequences)):
        new_hash = ""
        seq_ids.append(_sequences[i].id)
        while True:
            new_hash = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            if new_hash in hash_list:
                continue
            else:
                hash_list.append(new_hash)
                break
        _sequences[i].id = new_hash

    hash_map = []
    for i in range(len(hash_list)):
        hash_map.append((hash_list[i], seq_ids[i]))

    return [hash_map, _sequences]


def pull_recs(_arg1, _arg2):
    if os.path.isfile(_arg1):
        in_file = os.path.abspath(_arg1)
        substring = _arg2

    else:
        in_file = os.path.abspath(_arg2)
        substring = _arg1

    _output = []
    for _seq in SeqIO.parse(in_file, guess_format(in_file)):
        if _seq.description.find(substring) != -1:
            _output.append(_seq)

    return _output


def pull_seq_ends(_sequences, _amount, _which_end):
    seq_ends = []
    for _seq in _sequences:
        if _which_end == 'front':
            _seq.seq = _seq.seq[:_amount]

        elif _which_end == "rear":
            _seq.seq = _seq.seq[-1 * _amount:]

        else:
            sys.exit("Error: you much pick 'front' or 'rear' as the third argument in pull_seq_ends.")
        seq_ends.append(_seq)

    return seq_ends


def find_repeats(_sequences):
    unique_seqs = {}
    repeat_ids = {}
    repeat_seqs = {}

    # First find replicate IDs
    for _seq in _sequences:
        if _seq.id in repeat_ids:
            repeat_ids[_seq.id].append(_seq)
        elif _seq.id in unique_seqs:
            repeat_ids[_seq.id] = [_seq]
            repeat_ids[_seq.id].append(unique_seqs[_seq.id])
            del(unique_seqs[_seq.id])
        else:
            unique_seqs[_seq.id] = _seq

    # Then look for replicate sequences
    flip_uniqe = {}
    del_keys = []
    for key, value in unique_seqs.items():  # find and remove duplicates in/from the unique list
        value = str(value.seq)
        if value not in flip_uniqe:
            flip_uniqe[value] = [key]
        else:
            if value not in repeat_seqs:
                repeat_seqs[value] = [key]
                repeat_seqs[value] += flip_uniqe[value]
                if flip_uniqe[value][0] in unique_seqs:
                    del_keys.append(flip_uniqe[value][0])
            else:
                repeat_seqs[value].append(key)
            del_keys.append(unique_seqs[key])

    for key in del_keys:
        if key in unique_seqs:
            del(unique_seqs[key])

    for key, value in repeat_ids.items():  # find duplicates in the repeat ID list
        for blahh in value:
            blahh = str(blahh.seq)
            if blahh not in flip_uniqe:
                flip_uniqe[blahh] = [key]
            else:
                if blahh not in repeat_seqs:
                    repeat_seqs[blahh] = [key]
                    repeat_seqs[blahh] += flip_uniqe[blahh]

                else:
                    repeat_seqs[blahh].append(key)
    return [unique_seqs, repeat_ids, repeat_seqs]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="seq_tools.py", description="Commandline wrapper for all the fun functions in"
                                                                      "this file. Play with your sequences!")

    parser.add_argument('-ga', '--guess_alphabet', action='store')
    parser.add_argument('-gf', '--guess_format', action='store')
    parser.add_argument('-cs', '--clean_seq', action='store', help="")
    parser.add_argument('-tr', '--translate', action='store', help="Convert coding sequences into amino acid sequences")
    parser.add_argument('-ns', '--num_seqs', action='store',
                        help="Counts how many sequences are present in an input file")
    parser.add_argument('-cts', '--concat_seqs', action='store',
                        help="Concatenate a bunch of sequences into a single solid string.")
    parser.add_argument('-fd2p', '--map_features_dna2prot', action='store', nargs=2,
                        help="Arguments: <nucl_gb_file> <prot_file>")
    parser.add_argument('-fp2d', '--map_features_prot2dna', action='store', nargs=2,
                        help="Arguments: <prot_gb_file> <nucl_file>")
    parser.add_argument('-cf', '--combine_features', action='store', nargs=2, help="Arguments: <seq_file1> <seq_file2>")
    parser.add_argument('-cl', '--combine_files', action='store', nargs="+",
                        help="Arguments: <format> <files ... > <True|False (for mixing formats; optional)>")
    parser.add_argument('-sf', '--screw_formats', action='store', nargs=2,
                        help="Arguments: <in_file> <out_format>")
    parser.add_argument('-sfa', '--screw_formats_align', action='store', nargs=2,
                        help="Arguments: <in_file> <out_format>")
    parser.add_argument('-hsi', '--hash_seq_ids', action='store',
                        help="Rename all the identifiers in a sequence list to a 10 character hash.")
    parser.add_argument('-pr', '--pull_records', action='store', nargs=2,
                        help="Get all the records with ids containing a given string")
    parser.add_argument('-pe', '--pull_record_ends', action='store', nargs=3,
                        help="Get the ends (front or rear) of all sequences in a file."
                             "Arguments: <file> <amount (int)> <front|rear>")
    parser.add_argument('-fr', '--find_repeats', action='store',
                        help="Identify whether a file contains repeat sequences and/or sequence ids")

    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')
    parser.add_argument('-f', '--format', help="Some functions use this flag for output format", action='store')

    in_args = parser.parse_args()

    if in_args.format:
        out_format = in_args.format
    else:
        out_format = "fasta"

    # Translate CDS
    if in_args.translate:
        if guess_alphabet(in_args.translate) != "nucl":
            sys.exit("Error: you need to supply DNA or RNA sequences to translate")
        proteins = translate_cds(in_args.translate)

        for protein in proteins:
            print("%s\n" % protein.format(out_format))

    # Concatenate sequences
    if in_args.concat_seqs:
        print(concat_seqs(in_args.concat_seqs))

    # Count number of sequences in a file
    if in_args.num_seqs:
        print(count_sequences(in_args.num_seqs))

    # Find repeat sequences or ids
    if in_args.find_repeats:
        path = os.path.abspath(in_args.find_repeats)
        with open(path, "r") as ifile:
            sequences = list(SeqIO.parse(ifile, guess_format(path)))
        unique, rep_ids, rep_seqs = find_repeats(sequences)
        output = ""
        if len(rep_ids) > 0:
            output += "Records with duplicate IDs:\n"
            for next_id in rep_ids:
                output += "%s, " % next_id
            output = "%s\n\n" % output.strip(", ")

        else:
            output += "No records with duplicate IDs\n\n"
        if len(rep_seqs) > 0:
            output += "Records with duplicate sequences:\n"
            for next_id in rep_seqs:
                if len(rep_seqs[next_id]) > 1:
                    output += "("
                    for seq_id in rep_seqs[next_id]:
                        output += "%s, " % seq_id
                    output = "%s), " % output.strip(", ")
                else:
                    output += "%s, " % next_id
            output = "%s\n\n" % output.strip(", ")
        else:
            output += "No records with duplicate sequences\n\n"
        if len(unique) > 0:
            output += "Unique records:\n"
            for next_id in unique:
                output += "%s, " % next_id
            output = "%s" % output.strip(", ")
        else:
            output += "No unique records"
        print(output)

    # Pull sequence ends
    if in_args.pull_record_ends:
        path, amount, which_end = in_args.pull_record_ends
        amount = int(amount)
        with open(os.path.abspath(path), "r") as ifile:
            sequences = list(SeqIO.parse(ifile, guess_format(path)))
            new_seqs = pull_seq_ends(sequences, amount, which_end)
            for seq in new_seqs:
                seq.seq = _set_alphabet(seq.seq)
                print(seq.format(out_format))

    # Pull records
    if in_args.pull_records:
        arg1, arg2 = in_args.pull_records
        records = pull_recs(arg1, arg2)
        for rec in records:
            rec.seq = _set_alphabet(rec.seq)
            print(rec.format(out_format))

    # Hash sequence ids
    if in_args.hash_seq_ids:
        with open(os.path.abspath(in_args.hash_seq_ids), "r") as ifile:
            seq_format = guess_format(in_args.hash_seq_ids)
            seqs = list(SeqIO.parse(ifile, seq_format))
            hashed = hash_seqeunce_ids(seqs)
            print("# Hash table\n%s\n\n# Sequences\n" % hashed[0])
            for seq in hashed[1]:
                print(seq.format(seq_format))

    # Screw formats
    if in_args.screw_formats:
        with open(os.path.abspath(in_args.screw_formats[0]), "r") as ifile:
            seqs = SeqIO.parse(ifile, guess_format(in_args.screw_formats[0]))
            print(screw_formats(seqs, in_args.screw_formats[1]))

    # Screw formats align
    if in_args.screw_formats_align:
        with open(os.path.abspath(in_args.screw_formats_align[0]), "r") as ifile:
            align_format = guess_format(in_args.screw_formats_align[0])
            alignments = list(AlignIO.parse(ifile, align_format))

        print(screw_formats_align(alignments, in_args.screw_formats_align[1]))

    # Guess alphabet
    if in_args.guess_alphabet:
        print(guess_alphabet(in_args.guess_alphabet))

    # Clean Seq
    if in_args.clean_seq:
        seqs = clean_seq(in_args.clean_seq)
        output = ""
        for seq in seqs:
            output += "%s\n\n" % seq
        print(output.strip())

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