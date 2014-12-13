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
from random import sample, choice, randint
from math import ceil
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError


# ##################################################### WISH LIST #################################################### #
def get_genbank_file():
    x = 1


def run_prosite():
    x = 1


# ################################################ INTERNAL FUNCTIONS ################################################ #
def _set_alphabet(_sequences, alpha=None):  # update sequence alphabet in place
    if not alpha:
        alpha = guess_alphabet(_sequences)
    if alpha == "nucl":
        alpha = IUPAC.ambiguous_dna
    elif alpha == "prot":
        alpha = IUPAC.protein
    else:
        sys.exit("Error: Can't deterimine alphabet in _set_alphabet")
    for i in range(len(_sequences)):
        _sequences[i].seq.alphabet = alpha
    return _sequences


def _print_recs(_sequences):
    if len(_sequences) == 0:
        print("Nothing returned.", file=sys.stderr)
        return False
    _sequences = _set_alphabet(_sequences)

    if out_format == "phylipi":
        _output = phylipi(_sequences)

    elif out_format == "phylipis":
        _output = phylipi(_sequences, "strict")

    else:
        _output = ""
        for _rec in _sequences:
            try:
                _output += _rec.format(out_format) + "\n"
            except ValueError as e:
                print("Error: %s\n" % e, file=sys.stderr)

    if in_args.in_place and in_place_allowed:
        if not os.path.exists(in_args.sequence[0]):
            print("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a file. Nothing "
                  "was written.",
                  file=sys.stderr)
            print(_output.strip())
        else:
            with open(os.path.abspath(in_args.sequence[0]), "w") as ofile:
                ofile.write(_output)
            print("File over-written at:\n%s" % os.path.abspath(in_args.sequence[0]), file=sys.stderr)
    else:
        print(_output.strip())
# ################################################# HELPER FUNCTIONS ################################################# #


def sequence_list(sequence, _seq_format=None):  # Open a file and parse, or convert raw into a Seq object
    if isinstance(sequence, list):
        _sequences = sequence
    elif os.path.isfile(sequence):
        if not _seq_format:
            _seq_format = guess_format(sequence)
        if not _seq_format:
            sys.exit("Error: could not determine the format of your input sequence file. Explicitly set with -r flag.")
        with open(sequence, "r") as _infile:
            _sequences = list(SeqIO.parse(_infile, _seq_format))
    else:
        # dna_or_prot = IUPAC.protein if guess_alphabet(sequence) == "prot" else IUPAC.ambiguous_dna
        _sequences = [SeqRecord(Seq(sequence))]

    return _sequences


def phylipi(_sequences, _format="relaxed"):  # _format in ["strict", "relaxed"]
    max_id_length = 0
    max_seq_length = 0
    for _seq in _sequences:
        max_id_length = len(_seq.id) if len(_seq.id) > max_id_length else max_id_length
        max_seq_length = len(_seq.seq) if len(_seq.seq) > max_seq_length else max_seq_length

    _output = " %s %s\n" % (len(_sequences), max_seq_length)
    for _seq in _sequences:
        _seq_id = _seq.id.ljust(max_id_length) if _format == "relaxed" else _seq.id[:10].ljust(10)
        _output += "%s  %s\n" % (_seq_id, _seq.seq)

    return _output
# #################################################################################################################### #


def shuffle(_sequences):
    _sequences = sequence_list(_sequences)
    _output = []
    for _ in range(len(_sequences)):
        random_index = randint(1, len(_sequences)) - 1
        _output.append(_sequences.pop(random_index))
    return _output


def rna2dna(_sequences):
    _sequences = sequence_list(_sequences)
    _output = []
    for _seq in _sequences:
        _seq.seq = Seq(str(_seq.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
        _output.append(_seq)
    return _output


def dna2rna(_sequences):
    _sequences = sequence_list(_sequences)
    _output = []
    for _seq in _sequences:
        _seq.seq = Seq(str(_seq.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
        _output.append(_seq)
    return _output


def guess_alphabet(_sequences):  # Does not handle ambigious dna
    _sequences = sequence_list(_sequences)
    _sequence = ""
    for next_seq in _sequences:
        if len(_sequence) > 1000:
            break
        _sequence += re.sub("[NX\-?]", "", str(next_seq.seq))

    if len(_sequence) == 0:
        return None
    percent_dna = float(_sequence.count("A") + _sequence.count("G") + _sequence.count("T") +
                        _sequence.count("C") + _sequence.count("U")) / float(len(_sequence))
    if percent_dna > 0.95:
        return "nucl"
    else:
        return "prot"


def translate_cds(_sequences):
    _sequences = sequence_list(_sequences)
    _output = []
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
    _new_seq = ""
    concat_ids = []
    features = []
    alpha = None
    for _seq_list in _sequences:
        _seqs = sequence_list(_seq_list)
        if not alpha:
            alpha = guess_alphabet(_seqs)
        elif alpha != guess_alphabet(_seq_list):
            sys.exit("Error: You are trying to concatinate protein and nucleotide sequences together")

        for _seq in _seqs:
            location = FeatureLocation(len(_new_seq), len(_new_seq) + len(str(_seq.seq)))
            feature = SeqFeature(location=location, id=_seq.id, type=_seq.id[:15])
            features.append(feature)
            concat_ids.append(_seq.id)
            _new_seq += str(_seq.seq)

    concat_ids = "|".join(concat_ids)

    _output = [SeqRecord(Seq(_new_seq, alphabet=alpha), description=concat_ids, id="concatination", features=features)]
    return _output


def clean_seq(_sequences):  # from file or raw
    """remove fasta headers, numbers, and whitespace from sequence strings"""
    _sequences = sequence_list(_sequences)
    _output = []
    for _seq in _sequences:
        _seq.seq = re.sub(">.*", "", str(_seq.seq))
        _seq.seq = re.sub("[0-9\s\-\*]", "", str(_seq.seq))
        _seq.seq = str(_seq.seq).upper()
        _output.append(_seq)

    return _output  # returns a list of cleaned sequence objects


def guess_format(file_path):
    # currently just looks at file extension, but might want to get more fancy
    if not os.path.isfile(file_path):
        return None

    file_path = file_path.split(".")
    if file_path[-1] in ["fa", "fas", "fasta"]:
        return "fasta"
    if file_path[-1] in ["gb", "gp", "genbank"]:
        return "gb"
    if file_path[-1] in ["nex", "nxs", "nexus"]:
        return "nexus"
    if file_path[-1] in ["phy", "ph", "phylip"]:
        return "phylip-relaxed"
    if file_path[-1] in ["sto", "pfam", "stockholm"]:
        return "stockholm"
    return None


# Apply DNA features to protein sequences
def map_features_dna2prot(dna_seqs, prot_seqs):
    prot_dict = SeqIO.to_dict(prot_seqs)
    dna_dict = SeqIO.to_dict(dna_seqs)
    _new_seqs = {}
    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            print("Warning: %s is in protein file, but not cDNA file" % _seq_id, file=sys.stderr)
            continue

        _new_seqs[_seq_id] = prot_dict[_seq_id]

        for feature in dna_dict[_seq_id].features:
            start = (feature.location.start + 2) / 3
            end = feature.location.end / 3
            new_feature = SeqFeature(location=FeatureLocation(ceil(start), ceil(end)), type=feature.type)
            prot_dict[_seq_id].features.append(new_feature)

    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            print("Warning: %s is in cDNA file, but not protein file" % _seq_id, file=sys.stderr)

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    return _seqs_list


# Apply DNA features to protein sequences
def map_features_prot2dna(prot_seqs, dna_seqs):
    prot_dict = SeqIO.to_dict(prot_seqs)
    dna_dict = SeqIO.to_dict(dna_seqs)
    _new_seqs = {}
    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            print("Warning: %s is in protein file, but not cDNA file" % _seq_id, file=sys.stderr)
            continue

        _new_seqs[_seq_id] = dna_dict[_seq_id]

        for feature in prot_dict[_seq_id].features:
            start = feature.location.start * 3 - 2
            end = feature.location.end * 3
            new_feature = SeqFeature(location=FeatureLocation(start, end), type=feature.type)
            dna_dict[_seq_id].features.append(new_feature)

    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            print("Warning: %s is in cDNA file, but not protein file" % _seq_id, file=sys.stderr)

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    return _seqs_list


# Merge feature lists
def combine_features(seqs1, seqs2):  # These arguments are _sequence() lists
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
            sys.exit(error_mes)

    for _seq_id in seq_dict2:
        if type(seq_dict2[_seq_id].seq.alphabet) != type(reference_alphabet):
            error_mes = "You have mixed multiple alphabets into your sequences. Make sure everything is the same.\n" \
                        "\t%s in first set\n\tOffending alphabet: %s\n\tReference alphabet: %s" \
                        % (_seq_id, seq_dict2[_seq_id].seq.alphabet, reference_alphabet)
            sys.exit(error_mes)

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
        if re.search(_search, _seq.description) or re.search(_search, _seq.id) or re.search(_search, _seq.name):
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
            del_keys.append(unique_seqs[key].id)

    for key in del_keys:
        if key in unique_seqs:
            del(unique_seqs[key])

    for key, value in repeat_ids.items():  # find duplicates in the repeat ID list
        for _rep_seq in value:
            _rep_seq = str(_rep_seq.seq)
            if _rep_seq not in flip_uniqe:
                flip_uniqe[_rep_seq] = [key]
            else:
                if _rep_seq not in repeat_seqs:
                    repeat_seqs[_rep_seq] = [key]
                    repeat_seqs[_rep_seq] += flip_uniqe[_rep_seq]

                else:
                    repeat_seqs[_rep_seq].append(key)
    return [unique_seqs, repeat_ids, repeat_seqs]


def delete_records(_sequences, search_str):
    _new_seqs = []
    deleted = pull_recs(_sequences, search_str)
    for _seq in _sequences:
        if _seq in deleted:
            continue
        else:
            _new_seqs.append(_seq)
    return _new_seqs


def delete_repeats(_sequences, scope='all'):  # scope in ['all', 'ids', 'seqs']
    # First, remove duplicate IDs
    if scope in ['all', 'ids']:
        _unique, _rep_ids, _rep_seqs = find_repeats(_sequences)
        if len(_rep_ids) > 0:
            for _rep_id in _rep_ids:
                store_one_copy = pull_recs(_sequences, _rep_id)[0]
                _sequences = delete_records(_sequences, _rep_id)
                _sequences += [store_one_copy]

    # Then remove duplicate sequences
    if scope in ['all', 'seqs']:
        _unique, _rep_ids, _rep_seqs = find_repeats(_sequences)

        if len(_rep_seqs) > 0:
            _rep_seq_ids = []
            for _seq in _rep_seqs:
                _rep_seq_ids.append([])
                for _rep_seq_id in _rep_seqs[_seq]:
                    _rep_seq_ids[-1].append(_rep_seq_id)

            for _rep_seqs in _rep_seq_ids:
                for _rep_seq in _rep_seqs[1:]:
                    _sequences = delete_records(_sequences, _rep_seq)

    return _sequences


def rename(_sequences, query, replace=""):
    _sequences = sequence_list(_sequences)
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
    parser.add_argument('-fd2p', '--map_features_dna2prot', action='store_true',
                        help="Arguments: one cDNA file and one protein file")
    parser.add_argument('-fp2d', '--map_features_prot2dna', action='store_true',
                        help="Arguments: one cDNA file and one protein file")
    parser.add_argument('-ri', '--rename_ids', action='store', nargs=2,
                        help="Arguments: <pattern> <substitution>")
    parser.add_argument('-cf', '--combine_features', action='store_true',
                        help="Takes the features in two files and combines them for each sequence")
    parser.add_argument('-sf', '--screw_formats', action='store', help="Arguments: <out_format>")
    parser.add_argument('-sh', '--shuffle', action='store_true',
                        help="Randomly reorder the position of records in the file.")
    parser.add_argument('-hsi', '--hash_seq_ids', action='store_true',
                        help="Rename all the identifiers in a sequence list to a 10 character hash.")
    parser.add_argument('-pr', '--pull_records', action='store',
                        help="Get all the records with ids containing a given string")
    parser.add_argument('-pe', '--pull_record_ends', action='store', nargs=2,
                        help="Get the ends (front or rear) of all sequences in a file."
                             "Arguments: <amount (int)> <front|rear>")
    parser.add_argument('-dr', '--delete_records', action='store', nargs="+",
                        help="Remove reocrds from a file. The deleted IDs are sent to stderr.")
    parser.add_argument('-drp', '--delete_repeats', action='store_true',
                        help="Strip repeat records (ids and/or identical sequences")
    parser.add_argument('-fr', '--find_repeats', action='store_true',
                        help="Identify whether a file contains repeat sequences and/or sequence ids")
    parser.add_argument("-mg", "--merge", help="Group a bunch of seq files together", action="store_true")

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')
    parser.add_argument('-o', '--out_format', help="Some functions use this flag for output format", action='store')
    parser.add_argument('-f', '--in_format', help="If the file extension isn't sane, specify the format", action='store')
    
    in_args = parser.parse_args()

    if in_args.out_format:
        out_format = in_args.out_format
    else:
        out_format = guess_format(in_args.sequence[0])

    in_place_allowed = False
    seqs = sequence_list(in_args.sequence[0], in_args.in_format)

    # Shuffle
    if in_args.shuffle:
        in_place_allowed = True
        _print_recs(shuffle(seqs))

    # Delete repeats
    if in_args.delete_repeats:
        in_place_allowed = True
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        unique, rep_ids, rep_seqs = find_repeats(seqs)
        stderr_output = ""
        if len(rep_ids) > 0:
            stderr_output += "# Records with duplicate ids deleted (first instance retained)\n"
            counter = 1
            for seq in rep_ids:
                stderr_output += "%s\t" % seq
                if counter % columns == 0:
                    stderr_output = "%s\n" % stderr_output.strip()
                counter += 1
            stderr_output += "\n\n"

        rep_seq_ids = []
        for seq in rep_seqs:
            rep_seq_ids.append([])
            for rep_seq_id in rep_seqs[seq]:
                rep_seq_ids[-1].append(rep_seq_id)

        if len(rep_seq_ids) > 0:
            stderr_output += "# Records with duplicate sequence deleted (first instance retained)\n"
            counter = 1
            for rep_seqs in rep_seq_ids:
                for rep_seq in rep_seqs[1:]:
                    stderr_output += "%s\t" % rep_seq
                    if counter % columns == 0:
                        stderr_output = "%s\n" % stderr_output.strip()
                    counter += 1
            stderr_output += "\n"

        if stderr_output != "":
            print("# ################################################################ #", file=sys.stderr)
            print(stderr_output.strip(), file=sys.stderr)
            print("# ################################################################ #\n", file=sys.stderr)
            _print_recs(delete_repeats(seqs))

        else:
            print("No duplicate records found", file=sys.stderr)

    # Delete records
    if in_args.delete_records:
        in_place_allowed = True
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        new_list = list(seqs)
        deleted_seqs = []
        for next_pattern in in_args.delete_records:

            deleted_seqs += pull_recs(new_list, next_pattern)
            new_list = delete_records(new_list, next_pattern)

        if len(deleted_seqs) > 0:
                counter = 1
                output = "# ####################### Deleted records ######################## #\n"
                for seq in deleted_seqs:
                    output += "%s\t" % seq.id
                    if counter % columns == 0:
                        output = "%s\n" % output.strip()
                    counter += 1
                output = "%s\n# ################################################################ #\n" % output.strip()
                print(output, file=sys.stderr)

        if len(deleted_seqs) == 0:
            print("# ################################################################ #", file=sys.stderr)
            print("# No sequence identifiers match %s" % ", ".join(in_args.delete_records), file=sys.stderr)
            print("# ################################################################ #\n", file=sys.stderr)

        _print_recs(new_list)

    # Merge
    if in_args.merge:
        new_list = []
        for infile in in_args.sequence:
            new_list += sequence_list(infile)
        _print_recs(new_list)

    # Screw formats
    if in_args.screw_formats:
        in_place_allowed = True
        out_format = in_args.screw_formats
        if in_args.in_place:  # Need to change the file extension
            os.remove(in_args.sequence[0])
            in_args.sequence[0] = ".".join(os.path.abspath(in_args.sequence[0]).split(".")[:-1]) + "." + out_format
            open(in_args.sequence[0], "w").close()

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
        for seq in seqs:
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
        hashed = hash_seqeunce_ids(seqs)
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
        print(guess_format(seqs))

    # Map features from cDNA over to protein
    if in_args.map_features_dna2prot:
        in_place_allowed = True
        file1, file2 = in_args.sequence[:2]

        file1 = sequence_list(file1)
        file2 = sequence_list(file2)

        if guess_alphabet(file1) == guess_alphabet(file2):
            sys.exit("Error: You must provide one DNA file and one protein file")

        if guess_alphabet(file1) == "prot":
            prot = file1
            dna = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            prot = file2
            dna = file1

        new_seqs = map_features_dna2prot(dna, prot)
        out_format = "gb"
        _print_recs(new_seqs)

    # Map features from protein over to cDNA
    if in_args.map_features_prot2dna:
        in_place_allowed = True
        file1, file2 = in_args.sequence[:2]

        file1 = sequence_list(file1)
        file2 = sequence_list(file2)

        if guess_alphabet(file1) == guess_alphabet(file2):
            sys.exit("Error: You must provide one DNA file and one protein file")

        if guess_alphabet(file1) == "nucl":
            dna = file1
            prot = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            dna = file2
            prot = file1

        new_seqs = map_features_prot2dna(prot, dna)
        out_format = "gb"
        _print_recs(new_seqs)

    # Combine feature sets from two files into one
    if in_args.combine_features:
        file1, file2 = in_args.sequence[:2]
        file1 = sequence_list(file1)
        file2 = sequence_list(file2)
        new_seqs = combine_features(file1, file2)
        _print_recs(new_seqs)