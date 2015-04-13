#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 9 2015 

"""
Run all-by-all AlignMe
"""

import sys
import os
import re
import shutil
from MyFuncs import run_multicore_function, TempDir
from SeqBuddy import SeqBuddy, delete_metadata
from pssm import PSSM
import argparse
from multiprocessing import Lock
from subprocess import Popen
from Bio import SeqIO, AlignIO


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
    
    Popen("bloctopus %s/NameFile.txt %s %s %s %s %s %s -P > /dev/null 2>&1" %
          (_tmp_dir, _tmp_dir, out_dir, shutil.which("blastall"), shutil.which("blastpgp"), blastdb,
           shutil.which("makemat")), shell=True).wait()

    Popen("octopus %s/NameFile.txt %s/PSSM_PRF_FILES %s/RAW_PRF_FILES %s -N > /dev/null 2>&1" %
          (_tmp_dir, out_dir, out_dir, out_dir), shell=True).wait()

    # do a little reformatting of the .nnprf files --> could probably re-implement this in python...
    # I just grabbed this perl line from the web
    Popen("perl -e 'local $/; $_ = <>; s/END(.*)//gs; print' %s/%s.nnprf > %s/%s.temp; mv %s/%s.temp %s/%s.nnprf" %
          (nnprf_dir, seq_obj.id, nnprf_dir, seq_obj.id, nnprf_dir, seq_obj.id, nnprf_dir,
           seq_obj.id), shell=True).wait()

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

    clean_up(["%s.*" % output_loc])
    return


def score_alignme(alignme_file):
    with open(alignme_file, "r") as ifile:
        file_lines = ifile.readlines()

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

        # Also need to include some score from the PSSM...

    return round(tally / normalizing_len, 5) / 4  # The '4' is for the number of columns being compared


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
        base_name = "AlineMe_input"
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
    mafft_command = "einsi --thread -1 --quiet %s > %s/%s_einsi.fasta" % (seq_file, out_dir, base_name)
    print("\nMAFFT alignment\n%s" % mafft_command)
    Popen(mafft_command, shell=True).wait()
    print("  --> DONE")

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

    # Score all AlignMe runs, and push to file
    with open(alignme_scores_file, "w") as ofile:
        for pair in pairwise_array:
            try:
                score = score_alignme("%s/ALIGNME_FILES/%s-%s.prf" % (out_dir, pair[0], pair[1]))
                ofile.write("%s\t%s\t%s\n" % (pair[0], pair[1], score))
            except FileNotFoundError:
                pass