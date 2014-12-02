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
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError


# ################################################ INTERNAL FUNCTIONS ################################################ #
def _set_alphabet(_seq_list, alpha=None):  # update sequence alphabet in place
    if not alpha:
        alpha = guess_alphabet(_seq_list)
    if alpha == "nucl":
        alpha = IUPAC.ambiguous_dna
    elif alpha == "prot":
        alpha = IUPAC.protein
    else:
        sys.exit("Error: Can't deterimine alphabet in _set_alphabet")
    for i in range(len(_seq_list)):
        _seq_list[i].seq.alphabet = alpha
    return _seq_list


def _sequence_list(sequence):  # Open a file and parse, or convert raw into a Seq object
    if isinstance(sequence, list):
        _sequences = sequence
    elif os.path.isfile(sequence):
        _seq_format = guess_format(sequence)
        if not _seq_format:
            sys.exit("Error: could not determine the format of your input sequence file.")
        with open(sequence, "r") as infile:
            _sequences = list(SeqIO.parse(infile, _seq_format))
    else:
        # dna_or_prot = IUPAC.protein if guess_alphabet(sequence) == "prot" else IUPAC.ambiguous_dna
        _sequences = [SeqRecord(Seq(sequence))]

    return _sequences


def _print_recs(rec_list):
    rec_list = _set_alphabet(rec_list)
    _output = ""
    for _rec in rec_list:
        try:
            _output += _rec.format(out_format) + "\n"
        except ValueError as e:
            print("Error: %s\n" % e, file=sys.stderr)

    if in_args.in_place and in_place_allowed:
        if not os.path.exists(in_args.sequence[0]):
            print("Error: The -i flag was passed in, but the positional argument doesn't seem to be a file",
                  file=sys.stderr)
            print(_output.strip())
        else:
            with open(os.path.abspath(in_args.sequence[0]), "w") as ofile:
                ofile.write(_output)
            print("File over-written at:\n%s" % os.path.abspath(in_args.sequence[0]), file=sys.stderr)
    else:
        print(_output.strip())

# #################################################################################################################### #


def rna2dna(_sequences):
    _sequences = _sequence_list(_sequences)
    _output = []
    for _seq in _sequences:
        _seq.seq = Seq(str(_seq.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
        _output.append(_seq)
    return _output


def dna2rna(_sequences):
    _sequences = _sequence_list(_sequences)
    _output = []
    for _seq in _sequences:
        _seq.seq = Seq(str(_seq.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
        _output.append(_seq)
    return _output


def guess_alphabet(_sequences):  # Does not handle ambigious dna
    if not isinstance(_sequences, list):
        _sequences = _sequence_list(_sequences)

    _sequence = ""
    for next_seq in _sequences:
        if len(_sequence) > 1000:
            break
        _sequence += str(next_seq.seq)
    _sequence = re.sub("[NX]", "", _sequence)

    percent_dna = float(_sequence.count("A") + _sequence.count("G") + _sequence.count("T") +
                        _sequence.count("C") + _sequence.count("U")) / float(len(_sequence))
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


def concat_seqs(_sequences):  # TODO: Add each concatinated sequence as a record feature
    _output = ""
    concat_ids = ""
    for _seq_list in _sequences:
        _seq_list = _sequence_list(_seq_list)
        _sequences = [_seq.seq for _seq in clean_seq(_seq_list)]
        _id_list = [_seq.id for _seq in clean_seq(_seq_list)]
        _output += "".join(_sequences)
        concat_ids += "|".join(_id_list)
    alpha = guess_alphabet(_output)
    _output = [SeqRecord(Seq(_output, alphabet=alpha), description=concat_ids, id="concatination")]
    return _output


def clean_seq(_sequences):  # from file or raw
    """remove fasta headers, numbers, and whitespace from sequence strings"""
    _sequences = _sequence_list(_sequences)
    _output = []
    for _seq in _sequences:
        _seq.seq = re.sub(">.*", "", str(_seq.seq))
        _seq.seq = re.sub("[0-9\s\-\*]", "", str(_seq.seq))
        _seq.seq = str(_seq.seq).upper()
        _output.append(_seq)

    return _output  # returns a list of cleaned sequence objects


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
def combine_features(seqs1, seqs2):
    # make sure there are no repeat ids
    _unique, _rep_ids, _rep_seqs = find_repeats(seqs1)
    if len(_rep_ids) > 0:
        sys.exit("Error: There are repeat IDs in the first file provided\n%s" % _rep_ids)

    _unique, _rep_ids, _rep_seqs = find_repeats(seqs2)
    if len(_rep_ids) > 0:
        sys.exit("Error: There are repeat IDs in the second file provided\n%s" % _rep_ids)

    seq_dict1 = {}
    seq_dict2 = {}

    for _seq in seqs1:
        seq_dict1[_seq.id] = _seq

    for _seq in seqs2:
        seq_dict2[_seq.id] = _seq

    # make sure that we're comparing apples to apples across all sequences (i.e., same alphabet)
    reference_alphabet = sample(seq_dict1.items(), 1)[0][1].seq.alphabet
    for _seq_id in seq_dict1:
        if type(seq_dict1[_seq_id].seq.alphabet) != type(reference_alphabet):
            error_mes = "You have mixed multiple alphabets into your sequences. Make sure everything is the same.\n" \
                        "\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s" \
                        % (_seq_id, seq_dict1[_seq_id].seq.alphabet, reference_alphabet)
            print(error_mes, file=sys.stderr)
            return False

    for _seq_id in seq_dict2:
        if type(seq_dict2[_seq_id].seq.alphabet) != type(reference_alphabet):
            error_mes = "You have mixed multiple alphabets into your sequences. Make sure everything is the same.\n" \
                        "\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s" \
                        % (_seq_id, seq_dict2[_seq_id].seq.alphabet, reference_alphabet)
            print(error_mes, file=sys.stderr)
            return False

    _new_seqs = {}
    for _seq_id in seq_dict1:
        if _seq_id in seq_dict2:
            for feature in seq_dict2[_seq_id].features:
                seq_dict1[_seq_id].features.append(feature)
        else:
            print("Warning: %s is only in the first set of sequences" % _seq_id, file=sys.stderr)

        _new_seqs[_seq_id] = seq_dict1[_seq_id]

    for _seq_id in seq_dict2:
        if _seq_id not in seq_dict1:
            print("Warning: %s is only in the first set of sequences" % _seq_id, file=sys.stderr)
            _new_seqs[_seq_id] = seq_dict2[_seq_id]

    return [_new_seqs[_seq_id] for _seq_id in _new_seqs]


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


def pull_recs(_sequences, _search):
    _output = []
    for _seq in _sequences:
        if _seq.description.find(_search) != -1 or _seq.id.find(_search) != -1 or _seq.name.find(_search) != -1:
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


def rename(_sequences, query, replace=""):
    _sequences = _sequence_list(_sequences)
    _new_seqs = []
    for _seq in _sequences:
        new_name = re.sub(query, replace, _seq.id)
        _seq.id = new_name
        _seq.name = new_name
        _new_seqs.append(_seq)
    return _new_seqs


# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="seq_tools.py", description="Commandline wrapper for all the fun functions in"
                                                                      "this file. Play with your sequences!")

    parser.add_argument("sequence", help="Supply a file path or a raw sequence", nargs="+")

    parser.add_argument('-ga', '--guess_alphabet', action='store_true')
    parser.add_argument('-gf', '--guess_format', action='store_true')
    parser.add_argument('-cs', '--clean_seq', action='store_true', help="Strip out numbers and stuff")
    parser.add_argument('-tr', '--translate', action='store_true', help="Convert coding sequences into amino acid sequences")
    parser.add_argument('-d2r', '--transcribe', action='store_true', help="Convert DNA sequences to RNA")
    parser.add_argument('-r2d', '--back_transcribe', action='store_true', help="Convert RNA sequences to DNA")
    parser.add_argument('-li', '--list_ids', action='store_true',
                        help="Output all the sequence identifiers in a file. Use -p to specify # columns to write")
    parser.add_argument('-ns', '--num_seqs', action='store_true',
                        help="Counts how many sequences are present in an input file")
    parser.add_argument('-cts', '--concat_seqs', action='store_true',
                        help="Concatenate a bunch of sequences into a single solid string.")
    parser.add_argument('-fd2p', '--map_features_dna2prot', action='store', nargs=2,
                        help="Arguments: <nucl_gb_file> <prot_file>")  # Modify for new convention
    parser.add_argument('-fp2d', '--map_features_prot2dna', action='store', nargs=2,
                        help="Arguments: <prot_gb_file> <nucl_file>")  # Modify for new convention
    parser.add_argument('-ri', '--rename_ids', action='store', nargs=2,
                        help="Arguments: <pattern> <substitution>")
    parser.add_argument('-cf', '--combine_features', action='store_true',
                        help="Takes the features in two files and combines them for each sequence")
    parser.add_argument('-sf', '--screw_formats', action='store', help="Arguments: out_format>")
    parser.add_argument('-hsi', '--hash_seq_ids', action='store_true',
                        help="Rename all the identifiers in a sequence list to a 10 character hash.")
    parser.add_argument('-pr', '--pull_records', action='store',
                        help="Get all the records with ids containing a given string")
    parser.add_argument('-pe', '--pull_record_ends', action='store', nargs=2,
                        help="Get the ends (front or rear) of all sequences in a file."
                             "Arguments: <amount (int)> <front|rear>")
    parser.add_argument('-fr', '--find_repeats', action='store_true',
                        help="Identify whether a file contains repeat sequences and/or sequence ids")
    parser.add_argument("-mg", "--merge", help="Group a bunch of seq files together", action="store_true")

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')
    parser.add_argument('-f', '--format', help="Some functions use this flag for output format", action='store')
    
    in_args = parser.parse_args()

    if in_args.format:
        out_format = in_args.format
    else:
        out_format = "fasta"

    in_place_allowed = False
    seqs = _sequence_list(in_args.sequence[0])

    # Merge
    if in_args.merge:
        new_list = []
        for infile in in_args.sequence:
            new_list += _sequence_list(infile)
        _print_recs(new_list)

    # Screw formats
    if in_args.screw_formats:
        in_place_allowed = True
        out_format = in_args.screw_formats
        _print_recs(seqs)

    # Renaming
    if in_args.rename_ids:
        in_place_allowed = True
        seqs = rename(seqs, in_args.rename_ids[0], in_args.rename_ids[1])
        _print_recs(seqs)

    # Transcribe
    if in_args.transcribe:
        in_place_allowed = True
        if guess_alphabet(seqs) != "nucl":
            sys.exit("Error: You need to provide an unabmigious DNA sequence.")
        seqs = dna2rna(seqs)
        _print_recs(seqs)

    # Back Transcribe
    if in_args.back_transcribe:
        in_place_allowed = True
        if guess_alphabet(seqs) != "nucl":
            sys.exit("Error: You need to provide an unabmigious DNA sequence.")
        seqs = rna2dna(seqs)
        _print_recs(seqs)

    # List identifiers
    if in_args.list_ids:
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1
        output = ""
        counter = 1
        for seq in _sequence_list(seqs):
            output += "%s\t" % seq.id
            if counter % columns == 0:
                output = "%s\n" % output.strip()
            counter += 1
        print(output.strip())

    # Translate CDS
    if in_args.translate:
        in_place_allowed = True
        if guess_alphabet(seqs) != "nucl":
            sys.exit("Error: you need to supply DNA or RNA sequences to translate")
        _print_recs(translate_cds(seqs))

    # Concatenate sequences
    if in_args.concat_seqs:
        _print_recs(concat_seqs(in_args.sequence))

    # Count number of sequences in a file
    if in_args.num_seqs:
        print(len(seqs))

    # Find repeat sequences or ids
    if in_args.find_repeats:
        unique, rep_ids, rep_seqs = find_repeats(seqs)
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
        amount, which_end = in_args.pull_record_ends
        amount = int(amount)
        new_seqs = pull_seq_ends(seqs, amount, which_end)
        _print_recs(new_seqs)

    # Pull records
    if in_args.pull_records:
        search = in_args.pull_records
        records = pull_recs(seqs, search)
        _print_recs(records)

    # Hash sequence ids
    if in_args.hash_seq_ids:
        in_place_allowed = True
        hashed = hash_seqeunce_ids(_sequence_list(seqs))
        hash_table = "# Hash table\n"
        for seq in hashed[0]:
            hash_table += "%s,%s\n" % (seq[0], seq[1])
        print("%s\n" % hash_table, file=sys.stderr)
        _print_recs(hashed[1])

    # Guess alphabet
    if in_args.guess_alphabet:
        print(guess_alphabet(seqs))

    # Clean Seq
    if in_args.clean_seq:
        in_place_allowed = True
        seqs = clean_seq(seqs)
        output = ""
        for seq in seqs:
            output += "%s\n\n" % seq.seq
        print(output.strip())

    # Guess format
    if in_args.guess_format:
        print(guess_format(in_args.sequence[0]))

    # Map features from cDNA over to protein
    if in_args.map_features_dna2prot:
        in_place_allowed = True
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
        in_place_allowed = True
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
        file1, file2 = in_args.sequence[:2]
        file1 = _sequence_list(file1)
        file2 = _sequence_list(file2)
        new_seqs = combine_features(file1, file2)
        _print_recs(new_seqs)