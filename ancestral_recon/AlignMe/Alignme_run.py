#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 9 2015 

"""
Run all-by-all AlignMe

The program is a pipeline with 7 distinct steps, each of which can be the 'start' point of the run using the --resume/-r
parameter from the command line: 'octopus', 'top_file', 'mafft', 'pssm', 'psipred', 'alignme', or 'final_tally'.
Note that each step depends on the proper input from previous steps.
"""

import os
import re
import shutil
import argparse
from multiprocessing import Lock
from subprocess import Popen
from math import log
from MyFuncs import run_multicore_function, TempDir, DynamicPrint, SafetyValve, Timer, normalize
from SeqBuddy import SeqBuddy, delete_metadata, clean_seq
from pssm import PSSM
from Bio import SeqIO, AlignIO
from Bio.SubsMat import MatrixInfo, SeqMat


PHAT = SeqMat(MatrixInfo.phat75_73)
BLOSUM62 = SeqMat(MatrixInfo.blosum62)
ambiguous_X = {"A": 0, "R": -1, "N": -1, "D": -1, "C": -2, "Q": -1, "E": -1, "G": -1, "H": -1, "I": -1, "L": -1,
               "K": -1, "M": -1, "F": -1, "P": -2, "S": 0, "T": 0, "W": -2, "Y": -1, "V": -1}

for aa in ambiguous_X:
    pair = sorted((aa, "X"))
    pair = tuple(pair)
    PHAT[pair] = ambiguous_X[aa]
    BLOSUM62[pair] = ambiguous_X[aa]


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
    blastdb, _tmp_dir = args
    _tmp_dir = _tmp_dir.subdir()
    with open("%s/%s.fa" % (_tmp_dir, seq_obj.id), "w") as _ofile:
        SeqIO.write(seq_obj, _ofile, "fasta")

    with open(_tmp_dir + "/NameFile.txt", "w") as _ofile:
        seq_id = seq_obj.id.split(".")
        _ofile.write(seq_id[0] + "\n")

    # Sometimes octopus fails for some unknown reason. Try again if it does.
    valve = SafetyValve(5)
    while True:
        Popen("bloctopus %s/NameFile.txt %s %s %s %s %s %s -P > /dev/null 2>&1" %
              (_tmp_dir, _tmp_dir, out_dir, shutil.which("blastall"), shutil.which("blastpgp"), blastdb,
               shutil.which("makemat")), shell=True).wait()

        if os.path.isfile("%s/PSSM_PRF_FILES/%s.prf" % (out_dir, seq_obj.id)) and \
                os.path.isfile("%s/RAW_PRF_FILES/%s.prf" % (out_dir, seq_obj.id)):
            break
        else:
            valve.step("%s at bloctopus" % seq_obj.id)

    valve.global_reps = 5
    while True:
        Popen("octopus %s/NameFile.txt %s/PSSM_PRF_FILES %s/RAW_PRF_FILES %s -N > /dev/null 2>&1" %
              (_tmp_dir, out_dir, out_dir, out_dir), shell=True).wait()

        if os.path.isfile("%s/NN_PRF_FILES/%s.nnprf" % (out_dir, seq_obj.id)) and \
                os.path.isfile("%s/%s.top" % (out_dir, seq_obj.id)):
            break
        else:
            valve.step("%s at octopus" % seq_obj.id)

    # do a little reformatting of the .nnprf files --> Delete everything after the first group of values, before END 1
    with open("%s/%s.nnprf" % (nnprf_dir, seq_obj.id), "r") as _ifile:
        nnprf = _ifile.read()

    with open("%s/%s.nnprf" % (nnprf_dir, seq_obj.id), "w") as _ofile:
        regex = re.match(r'.+?(?=END)', nnprf, re.DOTALL)
        _ofile.write(regex.group(0))

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

    # AlignMe makes an error while writing alignment files if there is a dash (-) in the name, so repair the damage...
    with open("%s/%s-%s.aln" % (alignme_dir, combination[0], combination[1]), 'r') as _ifile:
        _output = _ifile.read()
        _output = re.sub("(^[A-Za-z]{3}\-.* )", r"\1 ", _output, flags=re.MULTILINE)

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


def mc_subs_mat_scores(_pair):
            align_file = "%s/ALIGNME_FILES/%s-%s.aln" % (out_dir, _pair[0], _pair[1])

            _align = SeqBuddy(align_file, "clustal")

            _align = _align.records
            subj_align = str(_align[0].seq)
            query_align = str(_align[1].seq)

            subj_top = str(top_seqs[_pair[0]].seq)
            query_top = str(top_seqs[_pair[1]].seq)

            output = alignment_sub_mat_score(subj_top, query_top, subj_align, query_align)
            with stdout_lock:
                with open(subs_mat_scores_file, "a") as _ofile:
                    _ofile.write("%s-%s\t%s\n" % (_pair[0], _pair[1], output["score"]))
            return


def mc_structural_scores(_pair):
            score = score_alignme("%s/ALIGNME_FILES/%s-%s.prf" % (out_dir, _pair[0], _pair[1]))
            with stdout_lock:
                with open(structural_scores_file, "a") as _ofile:
                    _ofile.write("%s-%s\t%s\n" % (_pair[0], _pair[1], score))

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
    parser.add_argument("-r", "--resume", help="Pick up at specific step a run left off, if possible.", choices=['octopus', 'top_file', 'mafft', 'pssm', 'psipred', 'alignme', 'final_tally'])

    in_args = parser.parse_args()

    temporary_directory = TempDir()
    timer = Timer()
    seqbuddy = []
    seq_set = ""

    print("\nPreparing SeqBuddy objects")
    for seq_set in in_args.input_sequences:
        seq_set = SeqBuddy(seq_set)
        seqbuddy += seq_set.records

    seqbuddy = SeqBuddy(seqbuddy)
    seqbuddy = delete_metadata(seqbuddy)
    seqbuddy = clean_seq(seqbuddy)

    expendable_dict = seqbuddy.to_dict()  # used during the creation of pairwise_array
    print("    ---> DONE")

    if not in_args.base_name:
        base_name = "Input"
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
    alignme_scores_file = "%s/%s_final_scores.csv" % (out_dir, base_name)
    top_file = "%s/%s_top.fasta" % (out_dir, base_name)
    seq_file = "%s/%s.fasta" % (out_dir, base_name)
    subs_mat_scores_file = "%s/%s_subsmat.csv" % (out_dir, base_name)
    structural_scores_file = "%s/%s_structural.csv" % (out_dir, base_name)

    if not in_args.resume:
        with open(seq_file, "w") as ofile:
            SeqIO.write(seqbuddy.records, ofile, "fasta")

        # unlink top file if it already exist...
        # These files are opened in append mode "a", so important to start with empty file
        clean_up([top_file])

    # multicore stuff
    stdout_lock = Lock()
    
    # Run octopus
    if not in_args.resume or in_args.resume == "octopus":
        print("\nExecuting Octopus")
        in_args.resume = False
        run_multicore_function(seqbuddy.records, octopus, ["/usr/local/blastdbs/pannexins",
                                                           temporary_directory])
        clean_up(["%s/PSSM_PRF_FILES" % out_dir, "%s/RAW_PRF_FILES" % out_dir, "%s/exec_times-bloctopus.txt" % out_dir,
                  "%s/exec_times-octopus.txt" % out_dir])

    # process top files
    if not in_args.resume or in_args.resume == "top_file":
        print("\nProcess top files")
        in_args.resume = False
        timer.start()
        top_file_handle = open(top_file, "a")
        for rec in seqbuddy.records:
            top_in_file = "%s/%s.top" % (out_dir, rec.id)
            with open(top_in_file, "r") as ifile:
                top_file_handle.write(ifile.read())
            clean_up([top_in_file])

        top_file_handle.close()
        print("    ---> DONE: %s" % timer.end())

    # create MSAs with MAFFT
    if not in_args.resume or in_args.resume == "mafft":
        print("\nMAFFT alignment running")
        in_args.resume = False
        timer.start()
        mafft_command = "einsi --thread -1 --quiet %s > %s/%s_einsi.fasta" % (seq_file, out_dir, base_name)
        Popen(mafft_command, shell=True).wait()
        print("    ---> DONE: %s" % timer.end())

    # set up all-by-all dict
    print("\nCreating all-by-all pairwise dictionary")
    timer.start()
    printer = DynamicPrint()
    dict_size = int(((len(seqbuddy.records) ** 2) - len(seqbuddy.records)) / 2)
    counter = 0
    pairwise_array = []
    for x in seqbuddy.records:
        del expendable_dict[x.id]
        for y in expendable_dict:
                pairwise_array.append((x.id, expendable_dict[y].id))
                counter += 1
                printer.write("%s of %s" % (counter, dict_size))

    print("\n    ---> DONE: %s" % timer.end())

    print("\nReading in top files")
    timer.start()
    with open(top_file, "r") as file:
        top_seqs = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    print("    ---> DONE: %s" % timer.end())

    print("\nReading in multiple sequence alignment")
    timer.start()
    with open("%s/%s_einsi.fasta" % (out_dir, base_name), "r") as seq_file:
        einsi_seqs = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    print("    ---> DONE: %s" % timer.end())

    # generate PSSMs -- differentiates between transmembrane and non-TM regions for
    # BLOSUM62 and PHAT respectively (implemented in the PSSM class).
    align = AlignIO.read("%s/%s_einsi.fasta" % (out_dir, base_name), "fasta")
    if not in_args.resume or in_args.resume == "pssm":
        in_args.resume = False
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
    if not in_args.resume or in_args.resume == "psipred":
        in_args.resume = False
        print("\nExecuting PSI-Pred")
        run_multicore_function(seqbuddy.records, _psi_pred)

    # Run AlignMe on all pairwise combinations
    if not in_args.resume or in_args.resume == "alignme":
        in_args.resume = False
        print("\nPairwise AlignMe runs:")
        run_multicore_function(pairwise_array, alignme)

    if not in_args.resume or in_args.resume == "final_tally":
        in_args.resume = False
        print("\nCompile substitution matrix scores.")
        # Get substitution matrix scores for each alignment
        with open(subs_mat_scores_file, "w"):
            pass

        run_multicore_function(pairwise_array, mc_subs_mat_scores)

        subs_mat_scores_handle = open(subs_mat_scores_file, "r")
        subs_mat_scores = {}
        for line in subs_mat_scores_handle:
            line = line.strip().split("\t")
            subs_mat_scores[line[0]] = float(line[1])

        subs_mat_scores_handle.close()
        subs_mat_scores = normalize(subs_mat_scores)

        # Score the structural components of all AlignMe runs
        print("\nCompile structural scores.")
        with open(structural_scores_file, "w"):
            pass

        run_multicore_function(pairwise_array, mc_structural_scores)
        structural_scores_handle = open(structural_scores_file, "r")
        structural_scores = {}

        for line in structural_scores_handle:
            line = line.strip().split("\t")
            structural_scores[line[0]] = float(line[1])

        structural_scores = normalize(structural_scores)

        # Combine scores and send to file
        print("\nPrinting final scores...")
        timer.start()
        with open(alignme_scores_file, "w") as ofile:
            for pair in pairwise_array:
                final_score = structural_scores["%s-%s" % (pair[0], pair[1])] + subs_mat_scores["%s-%s" % (pair[0], pair[1])]
                final_score /= 2
                ofile.write("%s\t%s\t%s\n" % (pair[0], pair[1], final_score))
        print("    ---> DONE: %s" % timer.end())