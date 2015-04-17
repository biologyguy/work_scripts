#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 9 2015 

"""
Run all-by-all AlignMe
"""

import os
import re
import shutil
import argparse
from multiprocessing import Lock
from subprocess import Popen
from math import log
from time import time
from MyFuncs import run_multicore_function, TempDir, pretty_time
from SeqBuddy import SeqBuddy, delete_metadata
from pssm import PSSM
from Bio import SeqIO, AlignIO
from Bio.SubsMat import MatrixInfo, SeqMat


PHAT = SeqMat(MatrixInfo.phat75_73)
BLOSUM62 = SeqMat(MatrixInfo.blosum62)


def clean_up(path_list):
    if type(path_list) != list:
        path_list = [path_list]
    for _next in path_list:
        if os.path.isdir(_next):
            # print("isdir: %s" % _next)
            shutil.rmtree(_next)
        elif os.path.isfile(_next):
            # print("isfile: %s" % _next)
            os.remove(_next)
    return


def octopus(seq_obj, args):
    valve = 1
    # There is a race condition (I think) that that causes a disruption in _tmp_dir.subdir(). It's infrequent
    # enough that I can just iterate until Octopus executes successfully.
    while True:
        blastdb, _tmp_dir = args
        try:
            _tmp_dir = _tmp_dir.subdir()
            with open("%s/%s.fa" % (_tmp_dir, seq_obj.id), "w") as _ofile:
                SeqIO.write(seq_obj, _ofile, "fasta")

            with open(_tmp_dir + "/NameFile.txt", "w") as _ofile:
                seq_id = seq_obj.id.split(".")
                _ofile.write(seq_id[0] + "\n")

            Popen("bloctopus %s/NameFile.txt %s %s %s %s %s %s -P > /dev/null 2>&1" %
                  (_tmp_dir, _tmp_dir, out_dir, shutil.which("blastall"), shutil.which("blastpgp"), blastdb,
                   shutil.which("makemat")), shell=True).wait()

            Popen("octopus %s/NameFile.txt %s/PSSM_PRF_FILES %s/RAW_PRF_FILES %s -N > /dev/null 2>&1" %
                  (_tmp_dir, out_dir, out_dir, out_dir), shell=True).wait()

            # do a little reformatting of the .nnprf files --> Delete everything after the first group of values, before END 1
            with open("%s/%s.nnprf" % (nnprf_dir, seq_obj.id), "r") as _ifile:
                nnprf = _ifile.read()

            with open("%s/%s.nnprf" % (nnprf_dir, seq_obj.id), "w") as _ofile:
                regex = re.match(r'.+?(?=END)', nnprf, re.DOTALL)
                _ofile.write(regex.group(0))
            break

        except FileNotFoundError:
            print("\nOops %s: %s" % (seq_obj.id, valve))
            valve += 1
            if valve == 20:
                raise SystemError("\nError: %s failed during Octopus step." % seq_obj.id)
            continue

    return


def _pssm(seq_obj, _pssm_lines):
    with open("%s/%s.pssm" % (pssm_dir, seq_obj.id), "w") as pssm_out_file:
        pssm_out_file.write("\nPSSM file for %s\n%s" % (seq_obj.id, _pssm_lines[0][1]))
        seq_list = list(seq_obj.seq)
        index = 2  # need to skip the first couple lines of the pssm file (header stuff)
        seq_position = 1
        for position in seq_list:
            if position != "-":
                out_line = re.search("[A-Z\-].+", _pssm_lines[0][index])
                pssm_out_file.write("%s %s%s\n" % (str(seq_position), position, out_line.group(0)[1:]))
                seq_position += 1
            index += 1
    return


def _psi_pred(seq_obj):
    tmp_fasta = "%s/%s.fa" % (tmp_dir, seq_obj.id)
    with open(tmp_fasta, "w") as tmp_file:
        SeqIO.write(seq_obj, tmp_file, "fasta")

    Popen("runpsipred %s > /dev/null 2>&1" % tmp_fasta, shell=True).wait()
    clean_up([seq_obj.id + ".ss", seq_obj.id + ".horiz"])
    Popen("mv %s.ss2 %s/%s.ss2" % (seq_obj.id, ss2_dir, seq_obj.id), shell=True).wait()
    return


def alignme(combination):
    sim_file_loc = "%s/%s-%s.simf" % (tmp_dir, combination[0], combination[1])

    # Magic numbser, might want to try running MCMCMC
    alignme_transmembrane_weight = 4.2
    alignme_psipred_weight = 1.4
    alignme_pssm_weight = 0.2

    with open(sim_file_loc, "w") as sim_file:
        sim_file.write("weight: %s type: UniversalProfileSimilarity column: 3 headerlines: 7 profile1: %s/%s.nnprf "
                       "profile2: %s/%s.nnprf\n" % (alignme_transmembrane_weight,
                                                    nnprf_dir, combination[0], nnprf_dir, combination[1]))

        ss2_profile_text = " headerlines: 2 profile1: %s/%s.ss2 profile2: %s/%s.ss2\n" % \
                           (ss2_dir, combination[0], ss2_dir, combination[1])
        sim_file.write("weight: %s type: UniversalProfileSimilarity column: 4%s" %
                       (alignme_psipred_weight, ss2_profile_text))
        sim_file.write("weight: %s type: UniversalProfileSimilarity column: 5%s" %
                       (alignme_psipred_weight, ss2_profile_text))
        sim_file.write("weight: %s type: UniversalProfileSimilarity column: 6%s" %
                       (alignme_psipred_weight, ss2_profile_text))

        sim_file.write("weight: %s type: PositionSpecificSimilarity PSSM1: %s/%s.pssm PSSM2: %s/%s.pssm\n" %
                       (alignme_pssm_weight, pssm_dir, combination[0], pssm_dir, combination[1]))

    output_loc = "%s/%s-%s" % (alignme_dir, combination[0], combination[1])

    # More magic numbers
    above_threshold_gap_opening_penalty = 5
    above_threshold_gap_extension_penalty = 3

    below_threshold_gap_opening_penalty = 3
    below_threshold_gap_extension_penalty = 1.5

    termini_gap_opening_penalty = 3
    termini_gap_extension_penalty = 1.5

    thresholds_for_penalties = 0.5

    strings = (tmp_dir, combination[0], tmp_dir, combination[1], sim_file_loc, output_loc, output_loc,
               below_threshold_gap_opening_penalty, above_threshold_gap_opening_penalty,
               below_threshold_gap_extension_penalty, above_threshold_gap_extension_penalty,
               termini_gap_opening_penalty, termini_gap_extension_penalty, thresholds_for_penalties)

    Popen("alignme -fasta_file1 %s/%s.fa -fasta_file2 %s/%s.fa -similarity_score_file %s -output_aligned_profiles "
          "%s.prf -output_aligned_sequences %s.aln -below_threshold_gap_opening_penalty %s "
          "-above_threshold_gap_opening_penalty %s -below_threshold_gap_extension_penalty %s "
          "-above_threshold_gap_extension_penalty %s -termini_gap_opening_penalty %s "
          "-termini_gap_extension_penalty %s -thresholds_for_penalties %s" % strings, shell=True).wait()

    # AlignMe makes an error while writing alignment files, so repair the damage...
    with open("%s/%s-%s.aln" % (alignme_dir, combination[0], combination[1]), 'r') as _ifile:
        _output = ''
        for line in _ifile:
            if line[0] == " ":
                line = line[1:]
            _output += line

    # AlignMe also writes an extra line to .aln files if the sequence's length is exactly divisible by line length. This
    # breaked BioPython...
    _output = re.sub("^[^ ].* +$", '', _output, flags=re.MULTILINE)

    with open("%s/%s-%s.aln" % (alignme_dir, combination[0], combination[1]), 'w') as _ofile:
        _ofile.write(_output)

    clean_up(["%s.*" % output_loc])
    return


def bit_score(raw_score):
    # These values were empirically determined for BLOSUM62 by Altschul
    bit_k_value = 0.035
    bit_lambda = 0.252

    bits = ((bit_lambda * raw_score) - (log(bit_k_value))) / log(2)
    return bits


def scale_range(data, percentile=1.0):
    if 0. > percentile > 1.0:
        raise ValueError("scale_range() percentile parameter should be between 0.5 and 1.0")

    if percentile < 0.5:
        percentile = 1 - percentile

    max_limit = round(len(data) * percentile) - 1
    min_limit = -1 * (max_limit + 1)

    if type(data) == dict:
        sorted_data = sorted([data[key] for key in data])
        _max = sorted_data[max_limit]
        _min = sorted_data[min_limit]
        data_range = _max - _min
        for key in data:
            data[key] = (data[key] - _min) / data_range
            data[key] = 1. if data[key] > 1. else data[key]
            data[key] = 0. if data[key] < 0. else data[key]

    else:
        sorted_data = sorted(data)
        _max = sorted_data[max_limit]
        _min = sorted_data[min_limit]
        data_range = _max - _min
        for i in range(len(data)):
            data[i] = (data[i] - _min) / data_range
            data[i] = 1. if data[i] > 1. else data[i]
            data[i] = 0. if data[i] < 0. else data[i]

    return data


def alignment_sub_mat_score(_subj_top, _query_top, _subj_align, _query_align):
    normalizing_len = (len(_subj_top) + len(_query_top)) / 2
    _score = 0
    gaps = 0
    for i in range(len(_subj_align)):
        if _subj_align[i] == "-":
            gaps += 1
            _query_top = _query_top[1:]

        elif _query_align[i] == "-":
            gaps += 1
            _subj_top = _subj_top[1:]

        else:
            _pair = sorted((_subj_align[i], _query_align[i]))
            _pair = tuple(_pair)
            if _subj_top[0] == "M" and _query_top[0] == "M":
                _score += PHAT[_pair]
            else:
                _score += BLOSUM62[_pair]

            _subj_top = _subj_top[1:]
            _query_top = _query_top[1:]

    _score = bit_score(_score)
    _score /= normalizing_len
    return {"score": _score, "gaps": gaps}


def score_alignme(alignme_file):  # This is the .prf file
    with open(alignme_file, "r") as _ifile:
        file_lines = _ifile.readlines()

    # clear out header rows
    while True:
        if not re.match("#", file_lines[0]):
            del file_lines[0]
        else:
            break
    del file_lines[-1]

    # determine normalizing maximum. I have semi-arbitrarily set this as one-half the sum of residues in both seqs
    seq1_len = 0.0
    seq2_len = 0.0
    for data in file_lines:
        data = re.split("\s+", data)
        if data[1] != "?0":
            seq1_len += 1
        if data[5] != "?0":
            seq2_len += 1

    normalizing_len = (seq1_len + seq2_len) / 2

    tally = 0.0
    for data in file_lines:
        if data[0] == "#":
            continue
        data = re.sub("\s+", ",", data)
        data = data.split(",")
        if "?0" in [data[1], data[5]]:  # Gaps get a score of 0
            continue

        # At the moment, all four parameters are given equal weight
        membrane_score = 1 - abs(float(data[1]) - float(data[5]))
        coil_score = 1 - abs(float(data[2]) - float(data[6]))
        helix_score = 1 - abs(float(data[3]) - float(data[7]))
        sheet_score = 1 - abs(float(data[4]) - float(data[8]))
        tally += membrane_score + coil_score + helix_score + sheet_score

    return round((tally / normalizing_len) / 4, 5)  # The '4' is for the number of columns being compared


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="Alignme_run", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="")

    parser.add_argument("input_sequences", nargs="+", action="store",
                        help="Select the sequence file(s) you want to analyze")
    parser.add_argument("-p", "--pssm", action="store",
                        help="Sepecify a PSSM file if pre-calculated. Otherwise it is calculate from the data.")
    parser.add_argument("-b", "--base_name", action="store", help="Specify the prefix for output.")
    parser.add_argument("-o", "--out_dir", default="./Alignme_output", action="store",
                        help="Where should output files go?")

    in_args = parser.parse_args()

    temporary_directory = TempDir()

    seqbuddy = []
    seq_set = ""

    for seq_set in in_args.input_sequences:
        seq_set = SeqBuddy(seq_set)
        seqbuddy += seq_set.records

    seqbuddy = SeqBuddy(seqbuddy)
    seqbuddy = delete_metadata(seqbuddy)

    if not in_args.base_name:
        base_name = "AlignMe_input"
    else:
        base_name = in_args.base_name

    # Output directories
    if not in_args.out_dir:
        out_dir = os.path.abspath("./Alignme_output")
    else:
        out_dir = os.path.abspath(in_args.out_dir)

    tmp_dir = temporary_directory.path
    pssm_dir = "%s/PSSM_FILES" % out_dir
    ss2_dir = "%s/SS2_FILES" % out_dir
    nnprf_dir = "%s/NN_PRF_FILES" % out_dir
    alignme_dir = "%s/ALIGNME_FILES" % out_dir
    
    # Create output directories
    for next_dir in [out_dir, pssm_dir, ss2_dir, nnprf_dir, alignme_dir]:
        if not os.path.exists(next_dir):
            os.makedirs(next_dir)
    
    # Output files
    alignme_scores_file = "%s/%s_scores.csv" % (out_dir, base_name)
    top_file = "%s/%s_top.fasta" % (out_dir, base_name)
    seq_file = "%s/%s.fasta" % (out_dir, base_name)

    with open(seq_file, "w") as ofile:
        SeqIO.write(seqbuddy.records, ofile, "fasta")

    # unlink top and trim files if they already exist...
    # These files are opened in append mode "a", so important to start with empty file
    clean_up([top_file])
    top_file_handle = open(top_file, "a")

    # multicore stuff
    stdout_lock = Lock()
    
    # Run octopus
    print("Execute Octopus on %s")
    run_multicore_function(seqbuddy.records, octopus, ["/usr/local/blastdbs/species_protein/Hydra_magnipapillata",
                                                       temporary_directory])
    clean_up(["%s/PSSM_PRF_FILES" % out_dir, "%s/RAW_PRF_FILES" % out_dir, "%s/exec_times-bloctopus.txt" % out_dir,
              "%s/exec_times-octopus.txt" % out_dir])

    # process top files
    print("\nProcess top files")
    for rec in seqbuddy.records:
        top_in_file = "%s/%s.top" % (out_dir, rec.id)
        with open(top_in_file, "r") as ifile:
            top_file_handle.write(ifile.read())
        clean_up([top_in_file])

    top_file_handle.close()
    print("  ---> DONE")

    # create MSAs with MAFFT
    mafft_time = round(time())
    mafft_command = "einsi --thread -1 --quiet %s > %s/%s_einsi.fasta" % (seq_file, out_dir, base_name)
    print("\nMAFFT alignment running")
    Popen(mafft_command, shell=True).wait()
    mafft_time = pretty_time(round(time()) - mafft_time)
    print("  --> DONE: %s" % mafft_time)

    # set up all-by-all dict
    pairwise_array = []
    for x in seqbuddy.records:
        for y in seqbuddy.records:
            if (x.id, y.id) in pairwise_array or (y.id, x.id) in pairwise_array or x.id == y.id:
                continue
            else:
                pairwise_array.append((x.id, y.id))
    
    with open(top_file, "r") as file:
        top_seqs = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    
    with open("%s/%s_einsi.fasta" % (out_dir, base_name), "r") as seq_file:
        einsi_seqs = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    
    # generate PSSMs -- differentiates between transmembrane and non-TM regions for
    # BLOSUM62 and PHAT respectively (implemented in the PSSM class).
    align = AlignIO.read("%s/%s_einsi.fasta" % (out_dir, base_name), "fasta")
    pssm = PSSM(align)
    pssm.name = base_name
    pssm.alignment_membranes(top_file)
    pssm.build_pssm()
    pssm.write("%s/%s.pssm" % (out_dir, base_name))

    # create custom pssms for each included sequence by deleting space rows
    print("\nCreating custom PSSMs for all sequences")
    with open("%s/%s.pssm" % (out_dir, base_name), "r") as pssm_file:
        pssm_lines = pssm_file.readlines()

    run_multicore_function(einsi_seqs, _pssm, [pssm_lines])

    # Now onto PSI-PRED
    print("\nExecuting PSI-Pred")
    run_multicore_function(seqbuddy.records, _psi_pred)

    # Run AlignMe on all pairwise combinations
    print("\nPairwise AlignMe runs:")
    run_multicore_function(pairwise_array, alignme)

    print("\nFinal scoring.")
    final_process_time = round(time())
    # Get substitution matrix scores for each alignment
    subs_mat_scores = {}
    for pair in pairwise_array:
        align = SeqBuddy("%s/ALIGNME_FILES/%s-%s.aln" % (out_dir, pair[0], pair[1]), "clustal")
        align = align.records
        subj_align = str(align[0].seq)
        query_align = str(align[1].seq)

        subj_top = str(top_seqs[pair[0]].seq)
        query_top = str(top_seqs[pair[1]].seq)

        output = alignment_sub_mat_score(subj_top, query_top, subj_align, query_align)
        subs_mat_scores["%s-%s" % (pair[0], pair[1])] = output["score"]
        # print("%s-%s\t%s\t%s" % (pair[0], pair[1], round(output["score"], 4), output["gaps"]))

    subs_mat_scores = scale_range(subs_mat_scores, 0.95)

    # Score the structural components of all AlignMe runs
    structural_scores = {}
    for pair in pairwise_array:
        try:
            score = score_alignme("%s/ALIGNME_FILES/%s-%s.prf" % (out_dir, pair[0], pair[1]))
            structural_scores["%s-%s" % (pair[0], pair[1])] = score

        except FileNotFoundError:
            print("Error: Failed to find %s/ALIGNME_FILES/%s-%s.prf" % (out_dir, pair[0], pair[1]))
            pass

    structural_scores = scale_range(structural_scores, 0.95)

    # Combine scores and send to file
    with open(alignme_scores_file, "w") as ofile:
        for pair in pairwise_array:
            final_score = structural_scores["%s-%s" % (pair[0], pair[1])] + subs_mat_scores["%s-%s" % (pair[0], pair[1])]
            final_score /= 2
            ofile.write("%s\t%s\t%s\n" % (pair[0], pair[1], final_score))

    final_process_time = pretty_time(round(time()) - final_process_time)
    print("  --> DONE: %s" % final_process_time)