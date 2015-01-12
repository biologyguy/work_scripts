#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 

"""
Collection of functions that do fun stuff with sequences. Pull them into a script, or run as a command line tool.
"""

# Standard library imports
import pdb
import sys
import os
import re
import string
from copy import copy
from random import sample, choice, randint
from math import ceil
from tempfile import TemporaryDirectory
from subprocess import Popen, PIPE
from shutil import which

# Third party package imports
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError
from scipy.cluster.hierarchy import linkage, to_tree, dendrogram
import numpy as np
import matplotlib.pylab as plt

# This will suppress the SearchIO warning, but be aware that new versions of BioPython may break SearchIO
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO


# ##################################################### WISH LIST #################################################### #
def get_genbank_file():
    x = 1
    return x


def run_prosite():
    x = 1
    return x


def order_features_by_alpha():
    x = 1
    return x


# ################################################# HELPER FUNCTIONS ################################################# #


class SeqBuddy():  # Open a file or read a handle and parse, or convert raw into a Seq object
    def __init__(self, _input, _in_format=None, _out_format=None):
        if not _in_format:
            self.in_format = guess_format(_input)
            self.out_format = str(self.in_format) if not _out_format else _out_format
        if not self.in_format:
            sys.exit("Error: could not determine the seq format in SeqBuddy(). "
                     "Try explicitly setting with -f flag.")

        self.out_format = self.in_format if not _out_format else _out_format

        if str(type(_input)) == "<class '__main__.SeqBuddy'>":
            _sequences = _input.seqs

        elif isinstance(_input, list):
            # make sure that the list is actually SeqIO records (just test a few...)
            for _seq in _input[:3]:
                if type(_seq) != SeqRecord:
                    sys.exit("Error: Seqlist is not populated with SeqRecords.")
            _sequences = _input

        elif str(type(_input)) == "<class '_io.TextIOWrapper'>":
            _sequences = list(SeqIO.parse(_input, self.in_format))

        elif os.path.isfile(_input):
            _input = open(_input, "r")
            _sequences = list(SeqIO.parse(_input, self.in_format))
            _input.close()
        else:
            _sequences = [SeqRecord(Seq(_input))]

        self.alpha = guess_alphabet(_sequences)

        for _i in range(len(_sequences)):
            _sequences[_i].seq.alphabet = self.alpha

        self.seqs = _sequences


def guess_alphabet(_seqs):  # Does not handle ambiguous dna
    _seqs = _seqs if isinstance(_seqs, list) else seqs.seqs
    _sequence = ""
    for next_seq in _seqs:
        if len(_sequence) > 1000:
            break
        _sequence += re.sub("[NX\-?]", "", str(next_seq.seq))
        _sequence = _sequence.upper()

    if len(_sequence) == 0:
        return None
    percent_dna = float(_sequence.count("A") + _sequence.count("G") + _sequence.count("T") +
                        _sequence.count("C") + _sequence.count("U")) / float(len(_sequence))
    if percent_dna > 0.95:
        nucl = IUPAC.ambiguous_rna if float(_sequence.count("U")) / float(len(_sequence)) > 0.05 else IUPAC.ambiguous_dna
        return nucl
    else:
        return IUPAC.protein


def guess_format(_input):  # _input can be list, SeqBuddy object, file handle, or file path.
    # If input is just a list, there is no BioPython in-format. Default to permissive gb.
    if isinstance(_input, list):
        return "gb"

    # Pull value directly from object if appropriate
    if str(type(_input)) == "<class '__main__.SeqBuddy'>":
        return _input.in_format

    # If input is a handle or path, try to read the file in each format, and assume success if not error and # seqs > 0
    if os.path.isfile(str(_input)):
        _input = open(_input, "r")

    if str(type(_input)) == "<class '_io.TextIOWrapper'>":
        possible_formats = ["phylip-relaxed", "stockholm", "fasta", "gb", "nexus"]
        for _format in possible_formats:
            try:
                _input.seek(0)
                _seqs = SeqIO.parse(_input, _format)
                if next(_seqs):
                    _input.seek(0)
                    return _format
                else:
                    continue
            except StopIteration:  # ToDo check that other types of error are not possible
                continue
            except ValueError:
                continue
        return None  # Unable to determine format from file handle

    else:
        sys.exit("Error: Unsupported _input argument in guess_format(). %s" % _input)


def phylipi(_input, _format="relaxed"):  # _format in ["strict", "relaxed"]
    max_id_length = 0
    max_seq_length = 0
    for _seq in _input.seqs:
        max_id_length = len(_seq.id) if len(_seq.id) > max_id_length else max_id_length
        max_seq_length = len(_seq.seq) if len(_seq.seq) > max_seq_length else max_seq_length

    _output = " %s %s\n" % (len(_input.seqs), max_seq_length)
    for _seq in _input.seqs:
        _seq_id = _seq.id.ljust(max_id_length) if _format == "relaxed" else _seq.id[:10].ljust(10)
        _output += "%s  %s\n" % (_seq_id, _seq.seq)

    return _output
# #################################################################################################################### #


def blast(_seqs, blast_db, blast_path=None, blastdbcmd=None):  # ToDo: Allow weird binary names to work
    if not blast_path:
        blast_path = which("blastp") if _seqs.alpha == IUPAC.protein else which("blastn")

    blast_check = Popen("%s -version" % blast_path, stdout=PIPE, shell=True).communicate()
    blast_check = re.search("([a-z])*[^:]", blast_check[0].decode("utf-8"))
    if blast_check:
        blast_check = blast_check.group(0)

    # ToDo Check NCBI++ tools are a conducive version (2.2.29 and above, I think [maybe .28])

    # Check to make sure blast is in path and ensure that the blast_db is present
    blast_db = os.path.abspath(blast_db)
    if blast_check == "blastp":
        if not which(blast_path):
            raise FileNotFoundError("blastp")

        if not os.path.isfile("%s.pin" % blast_db) or not os.path.isfile("%s.phr" % blast_db) \
                or not os.path.isfile("%s.psq" % blast_db):
            sys.exit("Error:\tBlastp database not found at '%s'" % blast_db)
    elif blast_check == "blastn":
        if not which(blast_path):
            raise FileNotFoundError("blastn")

        if not os.path.isfile("%s.nin" % blast_db) or not os.path.isfile("%s.nhr" % blast_db) \
                or not os.path.isfile("%s.nsq" % blast_db):
            sys.exit("Error:\tBlastn database not found at '%s'" % blast_db)
    else:
        sys.exit("Blast binary doesn't seem to work, at %s" % blast_path)

    if not blastdbcmd:
        blastdbcmd = "blastdbcmd"

    if not which(blastdbcmd):
        raise FileNotFoundError("blastdbcmd")

    # Check that blastdb was made with the -parse_seqids flag
    extensions = ["pog", "psd", "psi"] if blast_check == "blastp" else ["nog", "nsd", "nsi"]
    if not os.path.isfile("%s.%s" % (blast_db, extensions[0])) or not \
            os.path.isfile("%s.%s" % (blast_db, extensions[1])) or not \
            os.path.isfile("%s.%s" % (blast_db, extensions[2])):
        sys.exit("Error: Incorrect blastdb. When making the blast database, please use the -parse_seqids flag.")

    tmp_dir = TemporaryDirectory()
    with open("%s/tmp.fa" % tmp_dir.name, "w") as ofile:
        SeqIO.write(_seqs.seqs, ofile, "fasta")

    Popen("%s -db %s -query %s/tmp.fa -out %s/out.txt -num_threads 20 -evalue 0.01 -outfmt 6" %
          (blast_path, blast_db, tmp_dir.name, tmp_dir.name), shell=True).wait()

    with open("%s/out.txt" % tmp_dir.name, "r") as ifile:
        blast_results = SearchIO.parse(ifile, "blast-tab")
        _records = list(blast_results)

    hit_ids = []
    for record in _records:
        for hsp in record.hsps:
            hit_id = hsp.hit_id

            if hit_id in hit_ids:
                continue

            hit_ids.append(hit_id)

    ofile = open("%s/seqs.fa" % tmp_dir.name, "w")
    for hit_id in hit_ids:
        hit = Popen("blastdbcmd -db %s -entry 'lcl|%s'" % (blast_db, hit_id), stdout=PIPE, shell=True).communicate()
        hit = hit[0].decode("utf-8")
        hit = re.sub("lcl\|", "", hit)
        ofile.write("%s\n" % hit)

    ofile.close()

    with open("%s/seqs.fa" % tmp_dir.name, "r") as ifile:
        _new_seqs = SeqBuddy(ifile)

    return _new_seqs


def shuffle(_seqs):
    _output = []
    for _ in range(len(_seqs.seqs)):
        random_index = randint(1, len(_seqs.seqs)) - 1
        _output.append(_seqs.seqs.pop(random_index))
    _seqs.seqs = _output
    return _seqs


def rna2dna(_seqs):
    _output = []
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq.back_transcribe()), alphabet=IUPAC.ambiguous_dna)
        _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def dna2rna(_seqs):
    _output = []
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq.transcribe()), alphabet=IUPAC.ambiguous_rna)
        _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def translate_cds(_seqs):
    _output = []
    for _seq in _seqs.seqs:
        try:
            _seq.seq = _seq.seq.translate(cds=True, to_stop=True)
        except TranslationError as e1:
            _seq.seq = Seq(str(_seq.seq)[:(len(str(_seq.seq)) - len(str(_seq.seq)) % 3)])
            try:
                _seq.seq = _seq.seq.translate()
                sys.stderr.write("Warning: %s is not a standard CDS\t-->\t%s\n" % (_seq.id, e1))
            except TranslationError as e2:
                sys.stderr.write("Error: %s failed to translate\t-->\t%s\n" % (_seq.id, e2))

        _seq.seq.alphabet = IUPAC.protein
        _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def concat_seqs(_seqs):
    _new_seq = ""
    concat_ids = []
    features = []
    for _seq in _seqs.seqs:
        location = FeatureLocation(len(_new_seq), len(_new_seq) + len(str(_seq.seq)))
        feature = SeqFeature(location=location, id=_seq.id, type=_seq.id[:15])
        features.append(feature)
        concat_ids.append(_seq.id)
        _new_seq += str(_seq.seq)

    concat_ids = "|".join(concat_ids)
    _new_seq = [SeqRecord(Seq(_new_seq, alphabet=_seqs.alpha), description=concat_ids, id="concatination", features=features)]
    _seqs = SeqBuddy(_new_seq)
    _seqs.out_format = "gb"
    return _seqs


def clean_seq(_seqs):
    """remove all non-sequence chracters from sequence strings"""
    _output = []
    for _seq in _seqs.seqs:
        _seq.seq = str(_seq.seq).upper()
        if _seqs.alpha == IUPAC.protein:
            _seq.seq = Seq(re.sub("[^ACDEFGHIKLMNPQRSTVWXY]", "", str(_seq.seq)), alphabet=_seqs.alpha)
        else:
            _seq.seq = Seq(re.sub("[^ATGCXNU]", "", str(_seq.seq)), alphabet=_seqs.alpha)

        _output.append(_seq)

    _seqs.seqs = _output
    return _seqs


def delete_metadata(_seqs):
    _new_seqs = []
    for _seq in _seqs.seqs:
        _new_seqs.append(SeqRecord(Seq(str(_seq.seq), alphabet=_seqs.alpha), id=_seq.id, name='', description=''))
    _seqs.seqs = _new_seqs
    return _seqs


# Apply DNA features to protein sequences
def map_features_dna2prot(dna_seqs, prot_seqs):
    prot_dict = SeqIO.to_dict(prot_seqs.seqs)
    dna_dict = SeqIO.to_dict(dna_seqs.seqs)
    _new_seqs = {}
    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            sys.stderr.write("Warning: %s is in protein file, but not cDNA file\n" % _seq_id)
            continue

        _new_seqs[_seq_id] = prot_dict[_seq_id]

        for feature in dna_dict[_seq_id].features:
            start = feature.location.start / 3
            end = feature.location.end / 3
            location = FeatureLocation(ceil(start), ceil(end))
            feature.location = location
            prot_dict[_seq_id].features.append(feature)

    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            sys.stderr.write("Warning: %s is in cDNA file, but not protein file\n" % _seq_id)

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    _seqs = SeqBuddy(_seqs_list)
    _seqs.out_format = "gb"
    return _seqs


# Apply DNA features to protein sequences
def map_features_prot2dna(prot_seqs, dna_seqs):
    prot_dict = SeqIO.to_dict(prot_seqs.seqs)
    dna_dict = SeqIO.to_dict(dna_seqs.seqs)
    _new_seqs = {}
    for _seq_id in prot_dict:
        if _seq_id not in dna_dict:
            sys.stderr.write("Warning: %s is in protein file, but not cDNA file\n" % _seq_id)
            continue

        _new_seqs[_seq_id] = dna_dict[_seq_id]

        for feature in prot_dict[_seq_id].features:
            start = feature.location.start * 3
            end = feature.location.end * 3
            location = FeatureLocation(start, end)
            feature.location = location
            dna_dict[_seq_id].features.append(feature)

    for _seq_id in dna_dict:
        if _seq_id not in prot_dict:
            sys.stderr.write("Warning: %s is in cDNA file, but not protein file\n" % _seq_id)

    _seqs_list = [_new_seqs[_seq_id] for _seq_id in _new_seqs]
    _seqs = SeqBuddy(_seqs_list)
    _seqs.out_format = "gb"
    return _seqs


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

    for _seq in seqs1.seqs:
        seq_dict1[_seq.id] = _seq

    for _seq in seqs2.seqs:
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
            _seq_feats1 = []  # Test list so features common to both records are not duplicated
            for feature in seq_dict1[_seq_id].features:
                _seq_feats1.append("%s-%s-%s" % (feature.location.start, feature.location.end, feature.type))
            for feature in seq_dict2[_seq_id].features:
                feature_check = "%s-%s-%s" % (feature.location.start, feature.location.end, feature.type)
                if feature_check in _seq_feats1:
                    continue
                else:
                    seq_dict1[_seq_id].features.append(feature)
        else:
            sys.stderr.write("Warning: %s is only in the first set of sequences\n" % _seq_id)

        _new_seqs[_seq_id] = seq_dict1[_seq_id]

    for _seq_id in seq_dict2:
        if _seq_id not in seq_dict1:
            sys.stderr.write("Warning: %s is only in the first set of sequences\n" % _seq_id)
            _new_seqs[_seq_id] = seq_dict2[_seq_id]

    _new_seqs = SeqBuddy([_new_seqs[_seq_id] for _seq_id in _new_seqs], _out_format=seqs1.in_format)
    return _new_seqs


def order_features_by_position(_seqs):
    for _seq in _seqs.seqs:
        new_feature_list = []  # 2D list, with each item being [location, feature_obj]
        for _feature in _seq.features:
            feat_location = str(_feature.location)
            feat_location = int(re.search("([\d]+)", feat_location).group(0))
            new_feature = [feat_location, _feature]
            new_feature_index = 0

            for _i in range(len(new_feature_list)):
                if feat_location >= new_feature_list[_i][0]:
                    new_feature_index = _i
                    break

            new_feature_list.insert(new_feature_index, new_feature)

        _seq.features = [x[1] for x in new_feature_list]
        _seq.features.reverse()
    return _seqs


def hash_seqeunce_ids(_seqs):
    hash_list = []
    seq_ids = []
    for i in range(len(_seqs.seqs)):
        new_hash = ""
        seq_ids.append(_seqs.seqs[i].id)
        while True:
            new_hash = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            if new_hash in hash_list:
                continue
            else:
                hash_list.append(new_hash)
                break
        _seqs.seqs[i].id = new_hash
        _seqs.seqs[i].name = new_hash

    hash_map = []
    for i in range(len(hash_list)):
        hash_map.append((hash_list[i], seq_ids[i]))

    return [hash_map, _seqs]


def pull_recs(_seqs, _search):
    _output = []
    for _seq in _seqs.seqs:
        if re.search(_search, _seq.description) or re.search(_search, _seq.id) or re.search(_search, _seq.name):
            _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def pull_seq_ends(_seqs, _amount, _which_end):
    seq_ends = []
    for _seq in _seqs.seqs:
        if _which_end == 'front':
            _seq.seq = _seq.seq[:_amount]

        elif _which_end == "rear":
            _seq.seq = _seq.seq[-1 * _amount:]

        else:
            sys.exit("Error: you much pick 'front' or 'rear' as the third argument in pull_seq_ends.")
        seq_ends.append(_seq)
    _seqs.seqs = seq_ends
    return _seqs


def find_repeats(_seqs):
    unique_seqs = {}
    repeat_ids = {}
    repeat_seqs = {}

    # First find replicate IDs
    for _seq in _seqs.seqs:
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


def delete_records(_seqs, search_str):
    _output = []
    _deleted = pull_recs(copy(_seqs), search_str).seqs
    for _seq in _seqs.seqs:
        if _seq in _deleted:
            continue
        else:
            _output.append(_seq)
    _seqs.seqs = _output
    return _seqs


def delete_features(_seqs, _pattern):
    for _seq in _seqs.seqs:
        retained_features = []
        for _feature in _seq.features:
            if not re.search(_pattern, _feature.type):
                retained_features.append(_feature)
        _seq.features = retained_features
    return _seqs


def delete_repeats(_seqs, scope='all'):  # scope in ['all', 'ids', 'seqs']
    # First, remove duplicate IDs
    if scope in ['all', 'ids']:
        _unique, _rep_ids, _rep_seqs = find_repeats(_seqs)
        if len(_rep_ids) > 0:
            for _rep_id in _rep_ids:
                store_one_copy = pull_recs(copy(_seqs), _rep_id).seqs[0]
                delete_records(_seqs, _rep_id)
                _seqs.seqs.append(store_one_copy)

    # Then remove duplicate sequences
    if scope in ['all', 'seqs']:
        _unique, _rep_ids, _rep_seqs = find_repeats(_seqs)
        if len(_rep_seqs) > 0:
            _rep_seq_ids = []
            for _seq in _rep_seqs:
                _rep_seq_ids.append([])
                for _rep_seq_id in _rep_seqs[_seq]:
                    _rep_seq_ids[-1].append(_rep_seq_id)

            for _rep_seqs in _rep_seq_ids:
                for _rep_seq in _rep_seqs[1:]:
                    delete_records(_seqs, _rep_seq)

    return _seqs


def rename(_seqs, query, replace=""):  # TODO Allow a replacement pattern increment (like numbers)
    for _seq in _seqs.seqs:
        new_name = re.sub(query, replace, _seq.id)
        _seq.id = new_name
        _seq.name = new_name
    return _seqs


def purge(_seqs, threshold):  # ToDo: Implement a way to return a certain # of sequences (i.e. auto-determine threshold)
    purge_set = {}
    for _seq in _seqs.seqs:
        _unique = True
        for _seq_id in purge_set:
            purge_seq = purge_set[_seq_id]["seq"]
            blast_seqs = SeqBuddy([_seq, purge_seq])
            _blast_res = bl2seq(blast_seqs)
            bit_score = float(_blast_res.split("\t")[5])
            if bit_score >= threshold:
                _unique = False
                purge_set[_seq_id]["sim_seq"].append(_seq.id)

        if _unique:
            purge_set[_seq.id] = {"seq": _seq, "sim_seq": []}

    _output = {"seqs": [], "deleted": {}}
    for _seq_id in purge_set:
        _output["seqs"].append(purge_set[_seq_id]["seq"])
        _output["deleted"][_seq_id] = purge_set[_seq_id]["sim_seq"]

    _seqs.seqs = _output["seqs"]
    return [_seqs, _output["deleted"]]


def bl2seq(_seqs):  # Does an all-by-all analysis, and does not return sequences
    blast_bin = "blastp" if _seqs.alpha == IUPAC.protein else "blastn"
    if not which(blast_bin):
        sys.exit("Error: %s not present in $PATH.")  # ToDo: Implement -p flag

    tmp_dir = TemporaryDirectory()
    _seqs_copy = _seqs.seqs[1:]
    subject_file = "%s/subject.fa" % tmp_dir.name
    query_file = "%s/query.fa" % tmp_dir.name
    _output = ""
    for subject in _seqs.seqs:
        with open(subject_file, "w") as ifile:
            SeqIO.write(subject, ifile, "fasta")

        for query in _seqs_copy:
            with open(query_file, "w") as ifile:
                SeqIO.write(query, ifile, "fasta")

            _blast_res = Popen("%s -subject %s -query %s -outfmt 6" %
                               (blast_bin, subject_file, query_file), stdout=PIPE, shell=True).communicate()
            _blast_res = _blast_res[0].decode().split("\n")[0].split("\t")
            if len(_blast_res) == 1:
                _output += "%s\t%s\t0\t0\t0\t0\n" % (subject.id, query.id)
            else:
                # values are: query, subject, %_ident, length, evalue, bit_score
                _output += "%s\t%s\t%s\t%s\t%s\t%s\n" % (_blast_res[0], _blast_res[1], _blast_res[2],
                                                         _blast_res[3], _blast_res[10], _blast_res[11].strip())

        _seqs_copy = _seqs_copy[1:]
    return _output.strip()


def uppercase(_seqs):
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq).upper(), alphabet=_seq.seq.alphabet)
    return _seqs


def lowercase(_seqs):
    for _seq in _seqs.seqs:
        _seq.seq = Seq(str(_seq.seq).lower(), alphabet=_seq.seq.alphabet)
    return _seqs


def denroblast(_seqs):  # This does not work yet... See Kelly and Maini, 2013, PlosONE
    _blast_res = bl2seq(_seqs).split("\n")

    # format the data into a dictionary for easier manipulation
    dist_dict = {}
    for pair in _blast_res:
        pair = pair.split("\t")
        if pair[0] in dist_dict:
            dist_dict[pair[0]][pair[1]] = float(pair[5])
        else:
            dist_dict[pair[0]] = {pair[1]: float(pair[5])}

        if pair[1] in dist_dict:
            dist_dict[pair[1]][pair[0]] = float(pair[5])
        else:
            dist_dict[pair[1]] = {pair[0]: float(pair[5])}

    # make a numpy array and headings list
    dist_array = np.zeros([len(_seqs.seqs), len(_seqs.seqs)])
    headings = []
    i = 0
    for _seq1 in _seqs.seqs:
        headings.append(_seq1.id)
        j = 0
        for _seq2 in _seqs.seqs:
            if _seq1.id != _seq2.id:
                dist_array[i][j] = dist_dict[_seq1.id][_seq2.id]
            j += 1
        i += 1

    headings.reverse()
    print(headings)
    data_link = linkage(dist_array, method='complete')
    dendrogram(data_link, labels=headings)
    plt.xticks(fontsize=8, rotation=90)
    plt.savefig("dendrogram.svg", format='svg')
    plt.show()


# ################################################# COMMAND LINE UI ################################################## #
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="seq_tools.py", description="Commandline wrapper for all the fun functions in"
                                                                      " this file. Play with your sequences!")

    parser.add_argument("sequence", help="Supply a file path or a raw sequence", nargs="+", default=sys.stdin)

    parser.add_argument('-ga', '--guess_alphabet', action='store_true')
    parser.add_argument('-gf', '--guess_format', action='store_true')
    parser.add_argument('-cs', '--clean_seq', action='store_true',
                        help="Strip out non-sequence characters, such as stops (*) and gaps (-)")
    parser.add_argument('-uc', '--uppercase', action='store_true',
                        help='Convert all sequences to uppercase')
    parser.add_argument('-lc', '--lowercase', action='store_true',
                        help='Convert all sequences to lowercase')
    parser.add_argument('-dm', '--delete_metadata', action='store_true',
                        help="Remove meta-data from file (only id is retained)")
    parser.add_argument('-rs', '--raw_seq', action='store_true',
                        help="Return line break separated sequences")
    parser.add_argument('-tr', '--translate', action='store_true',
                        help="Convert coding sequences into amino acid sequences")
    parser.add_argument('-d2r', '--transcribe', action='store_true',
                        help="Convert DNA sequences to RNA")
    parser.add_argument('-r2d', '--back_transcribe', action='store_true',
                        help="Convert RNA sequences to DNA")
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
    parser.add_argument('-ri', '--rename_ids', action='store', metavar=('<pattern>', '<substitution>'), nargs=2,
                        help="Replace some pattern in ids with something else.")
    parser.add_argument('-cf', '--combine_features', action='store_true',
                        help="Takes the features in two files and combines them for each sequence")
    parser.add_argument('-ofp', '--order_features_by_position', action='store_true',
                        help="Change the output order of sequence features, based on sequence position")
    parser.add_argument('-sf', '--screw_formats', action='store', metavar="<out_format>",
                        help="Change the file format to something else.")
    parser.add_argument('-sh', '--shuffle', action='store_true',
                        help="Randomly reorder the position of records in the file.")
    parser.add_argument('-hsi', '--hash_seq_ids', action='store_true',
                        help="Rename all the identifiers in a sequence list to a 10 character hash.")
    parser.add_argument('-pr', '--pull_records', action='store',
                        help="Get all the records with ids containing a given string")
    parser.add_argument('-pe', '--pull_record_ends', action='store', nargs=2, metavar="<amount (int)> <front|rear>",
                        help="Get the ends (front or rear) of all sequences in a file.")
    parser.add_argument('-dr', '--delete_records', action='store', nargs="+",
                        help="Remove reocrds from a file. The deleted IDs are sent to stderr.")
    parser.add_argument('-df', '--delete_features', action='store', nargs="+",
                        help="Remove specified features from all records.")
    parser.add_argument('-drp', '--delete_repeats', action='store_true',
                        help="Strip repeat records (ids and/or identical sequences")
    parser.add_argument('-fr', '--find_repeats', action='store_true',
                        help="Identify whether a file contains repeat sequences and/or sequence ids")
    parser.add_argument("-mg", "--merge", action="store_true",
                        help="Group a bunch of seq files together",)
    parser.add_argument("-bl", "--blast", metavar="BLAST database", action="store",
                        help="BLAST your sequence file using common settings, return the hits from blastdb")
    parser.add_argument("-bl2s", "--bl2seq", action="store_true",
                        help="All-by-all blast amoung sequences using bl2seq. Returns only the top hit from each search.")
    parser.add_argument("-prg", "--purge", metavar="Max BLAST score", type=int, action="store",
                        help="Delete sequences with high similarity")
    parser.add_argument("-drb", "--dendroblast", action="store_true",
                        help="Create a dendrogram from pairwise blast bit-scores. Returns newick format.")

    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')
    parser.add_argument('-p', '--params', help="Free form arguments for some functions", nargs="+", action='store')
    parser.add_argument('-o', '--out_format', help="If you want a specific format output", action='store')
    parser.add_argument('-f', '--in_format', help="If SeqBuddy can't guess the file format, just specify it directly.",
                        action='store')
    
    in_args = parser.parse_args()

    in_place_allowed = False

    seqs = []
    seq_set = ""
    for seq_set in in_args.sequence:
        seq_set = SeqBuddy(seq_set, in_args.in_format)
        seqs += seq_set.seqs

    seqs = SeqBuddy(seqs)

    seqs.out_format = in_args.out_format if in_args.out_format else seq_set.out_format

    # ############################################# INTERNAL FUNCTION ################################################ #
    def _print_recs(_seqs):
        if len(_seqs.seqs) == 0:
            sys.stderr.write("Nothing returned.\n")
            return False

        if _seqs.out_format == "phylipi":
            _output = phylipi(_seqs)

        elif _seqs.out_format == "phylipis":
            _output = phylipi(_seqs, "strict")

        else:
            tmp_dir = TemporaryDirectory()
            with open("%s/seqs.tmp" % tmp_dir.name, "w") as ofile:
                SeqIO.write(_seqs.seqs, ofile, _seqs.out_format)

            with open("%s/seqs.tmp" % tmp_dir.name, "r") as ifile:
                _output = ifile.read()

        if in_args.in_place and in_place_allowed:
            if not os.path.exists(in_args.sequence[0]):
                sys.stderr.write("Warning: The -i flag was passed in, but the positional argument doesn't seem to be a "
                                 "file. Nothing was written.\n")
                sys.stdout.write("%s\n" % _output.strip())
            else:
                with open(os.path.abspath(in_args.sequence[0]), "w") as ofile:
                    ofile.write(_output)
                sys.stderr.write("File over-written at:\n%s\n" % os.path.abspath(in_args.sequence[0]))
        else:
            sys.stdout.write("%s\n" % _output.strip())

    def _get_blast_binaries():
        blastp = None
        blastn = None
        blastdbcmd = None
        if in_args.params:
            for param in in_args.params:
                binary = Popen("%s -version" % param, stdout=PIPE, shell=True).communicate()
                binary = re.search("([a-z])*[^:]", binary[0].decode("utf-8"))
                binary = binary.group(0)
                if binary == "blastp":
                    blastp = param
                elif binary == "blastn":
                    blastn = param
                elif binary == "blastdbcmd":
                    blastdbcmd = param

        blastp = blastp if blastp else which("blastp")
        blastn = blastn if blastn else which("blastn")
        blastdbcmd = blastdbcmd if blastdbcmd else which("blastdbcmd")

        return {"blastdbcmd": blastdbcmd, "blastp": blastp, "blastn": blastn}

    # ############################################## COMMAND LINE LOGIC ############################################## #
    # DendroBlast
    if in_args.dendroblast:
        denroblast(seqs)

    # Purge
    if in_args.purge:
        in_place_allowed = True
        purged_seqs, deleted = purge(seqs, in_args.purge)

        stderr_output = "# Deleted record mapping #\n"
        for seq_id in deleted:
            stderr_output += "%s\n%s\n\n" % (seq_id, deleted[seq_id])

        sys.stderr.write(stderr_output)
        _print_recs(purged_seqs)

    # BL2SEQ
    if in_args.bl2seq:
        sys.stderr.write("#query\tsubject\t%_ident\tlength\tevalue\tbit_score\n")
        sys.stdout.write("%s\n" % bl2seq(seqs))

    # BLAST
    if in_args.blast:
        blast_binaries = _get_blast_binaries()
        blast_binary_path = blast_binaries["blastp"] if seqs.alpha == IUPAC.protein else blast_binaries["blastn"]
        try:
            blast_res = blast(seqs, in_args.blast, blast_path=blast_binary_path, blastdbcmd=blast_binaries["blastdbcmd"])

        except FileNotFoundError as e:
            sys.exit("%s binary not found, explicitly set with the -p flag.\n"
                     "To pass in the path to both blast(p/n) and blastdbcmd, separate them with a space." % e)
        _print_recs(blast_res)

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
            sys.stderr.write("# ################################################################ #\n")
            sys.stderr.write("%s\n" % stderr_output.strip())
            sys.stderr.write("# ################################################################ #\n")

            _print_recs(delete_repeats(seqs))

        else:
            sys.stderr.write("No duplicate records found\n")

    # Delete records
    if in_args.delete_records:
        in_place_allowed = True
        if in_args.params:
            columns = int(in_args.params[0])
        else:
            columns = 1

        new_list = SeqBuddy(list(seqs.seqs))
        deleted_seqs = []
        for next_pattern in in_args.delete_records:
            deleted_seqs += pull_recs(copy(new_list), next_pattern).seqs
            delete_records(new_list, next_pattern)

        if len(deleted_seqs) > 0:
                counter = 1
                output = "# ####################### Deleted records ######################## #\n"
                for seq in deleted_seqs:
                    output += "%s\t" % seq.id
                    if counter % columns == 0:
                        output = "%s\n" % output.strip()
                    counter += 1
                output = "%s\n# ################################################################ #\n" % output.strip()
                sys.stderr.write(output)

        if len(deleted_seqs) == 0:
            sys.stderr.write("# ################################################################ #\n")
            sys.stderr.write("# No sequence identifiers match %s\n" % ", ".join(in_args.delete_records))
            sys.stderr.write("# ################################################################ #\n")

        new_list.out_format = in_args.out_format if in_args.out_format else seqs.out_format
        _print_recs(new_list)

    # Delete features
    if in_args.delete_features:
        in_place_allowed = True
        for next_pattern in in_args.delete_features:
            delete_features(seqs, next_pattern)
        _print_recs(seqs)

    # Merge
    if in_args.merge:
        new_list = SeqBuddy([])
        for infile in in_args.sequence:
            new_list.seqs += SeqBuddy(infile).seqs

        new_list.out_format = in_args.out_format if in_args.out_format else seqs.out_format
        _print_recs(new_list)

    # Screw formats
    if in_args.screw_formats:
        in_place_allowed = True
        seqs.out_format = in_args.screw_formats
        if in_args.in_place:  # Need to change the file extension
            os.remove(in_args.sequence[0])
            in_args.sequence[0] = ".".join(os.path.abspath(in_args.sequence[0]).split(".")[:-1]) + "." + seqs.out_format
            open(in_args.sequence[0], "w").close()

        _print_recs(seqs)

    # Renaming
    if in_args.rename_ids:
        in_place_allowed = True
        seqs = rename(seqs, in_args.rename_ids[0], in_args.rename_ids[1])
        _print_recs(seqs)

    # Uppercase
    if in_args.uppercase:
        in_place_allowed = True
        _print_recs(uppercase(seqs))

    # Lowercase
    if in_args.lowercase:
        in_place_allowed = True
        _print_recs(lowercase(seqs))

    # Transcribe
    if in_args.transcribe:
        in_place_allowed = True
        if seqs.alpha != IUPAC.ambiguous_dna:
            sys.exit("Error: You need to provide a DNA sequence.")
        seqs = dna2rna(seqs)
        _print_recs(seqs)

    # Back Transcribe
    if in_args.back_transcribe:
        in_place_allowed = True
        if seqs.alpha != IUPAC.ambiguous_rna:
            sys.exit("Error: You need to provide an RNA sequence.")
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
        for seq in seqs.seqs:
            output += "%s\t" % seq.id
            if counter % columns == 0:
                output = "%s\n" % output.strip()
            counter += 1
        sys.stdout.write("%s\n" % output.strip())

    # Translate CDS
    if in_args.translate:
        in_place_allowed = True
        if seqs.alpha == IUPAC.protein:
            sys.exit("Error: you need to supply DNA or RNA sequences to translate")
        _print_recs(translate_cds(seqs))

    # Concatenate sequences
    if in_args.concat_seqs:
        seqs = concat_seqs(seqs)
        if in_args.out_format:
            seqs.out_format = in_args.out_format
        _print_recs(seqs)

    # Count number of sequences in a file
    if in_args.num_seqs:
        sys.stdout.write("%s\n" % len(seqs.seqs))

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
        sys.stdout.write("%s\n" % output)

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
        sys.stderr.write("%s\n" % hash_table)
        _print_recs(hashed[1])

    # Guess alphabet
    if in_args.guess_alphabet:
        if seqs.alpha == IUPAC.protein:
            sys.stdout.write("prot\n")
        elif seqs.alpha == IUPAC.ambiguous_dna:
            sys.stdout.write("dna\n")
        elif seqs.alpha == IUPAC.ambiguous_rna:
            sys.stdout.write("rna\n")
        else:
            sys.stdout.write("Undetermined\n")

    # Delete metadata
    if in_args.delete_metadata:
        in_place_allowed = True
        _print_recs(delete_metadata(seqs))

    # Raw Seq
    if in_args.raw_seq:
        seqs = clean_seq(seqs)
        output = ""
        for seq in seqs.seqs:
            output += "%s\n\n" % seq.seq
        sys.stdout.write("%s\n" % output.strip())

    # Clean Seq
    if in_args.clean_seq:
        in_place_allowed = True
        _print_recs(clean_seq(seqs))

    # Guess format
    if in_args.guess_format:
        for seq_set in in_args.sequence:
            sys.stdout.write("%s\t-->\t%s\n" % (seq_set, SeqBuddy(seq_set).in_format))

    # Map features from cDNA over to protein
    if in_args.map_features_dna2prot:
        in_place_allowed = True
        file1, file2 = in_args.sequence[:2]

        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)

        if file1.alpha == file2.alpha:
            sys.exit("Error: You must provide one DNA file and one protein file")

        if file1.alpha == IUPAC.protein:
            prot = file1
            dna = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            prot = file2
            dna = file1

        new_seqs = map_features_dna2prot(dna, prot)
        _print_recs(new_seqs)

    # Map features from protein over to cDNA
    if in_args.map_features_prot2dna:
        in_place_allowed = True
        file1, file2 = in_args.sequence[:2]

        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)

        if file1.alpha == file2.alpha:
            sys.exit("Error: You must provide one DNA file and one protein file")

        if file1.alpha != IUPAC.protein:
            dna = file1
            prot = file2
        else:
            in_args.sequence[0] = in_args.sequence[1]  # in case the -i flag is thrown
            dna = file2
            prot = file1

        new_seqs = map_features_prot2dna(prot, dna)
        _print_recs(new_seqs)

    # Combine feature sets from two files into one
    if in_args.combine_features:
        file1, file2 = in_args.sequence[:2]
        file1 = SeqBuddy(file1)
        file2 = SeqBuddy(file2)
        new_seqs = combine_features(file1, file2)
        _print_recs(new_seqs)

    if in_args.order_features_by_position:
        in_place_allowed = True
        _print_recs(order_features_by_position(seqs))