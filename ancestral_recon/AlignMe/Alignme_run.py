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
from SeqBuddy import SeqBuddy
from pssm import PSSM
import argparse
from multiprocessing import Lock
from subprocess import Popen
from Bio import SeqIO


def clean_up(path_list):
    if type(path_list) != list:
        path_list = [path_list]
    for _next in path_list:
        if os.path.isdir(_next):
            print("isdir: %s" % _next)
            shutil.rmtree(_next)
        elif os.path.isfile(_next):
            print("isfile: %s" % _next)
            os.remove(_next)
    return


def octopus(seq_obj, args):
    blastdb, _tmp_dir = args
    _tmp_dir = _tmp_dir.subdir()
    with open("%s/%s.fa" % (_tmp_dir, seq_obj.id), "w") as file:
        SeqIO.write(seq_obj, file, "fasta")

    with open(_tmp_dir + "/NameFile.txt", "w") as file:
        seq_id = seq_obj.id.split(".")
        file.write(seq_id[0] + "\n")
    
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="Alignme_run", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input_file", help="Select the sequence file you want to analyze", action="store")
    parser.add_argument("-p", "--pssm", action="store",
                        help="Sepecify a PSSM file if pre-calculated. Otherwise it is calculate from the data.")
    parser.add_argument("-b", "--base_name", action="store", help="Specify the prefix for output.")
    parser.add_argument("-o", "--out_dir", default="./Alignme_output", action="store", help="Where should output files go?")

    in_args = parser.parse_args()

    temporary_directory = TempDir()
    seqbuddy = SeqBuddy(in_args.input_file)

    if not in_args.base_name:
        base_name = ".".join(in_args.input_file.split("/")[-1].split(".")[:-1])
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
    top_file = "%s/%s_top.fasta" % (out_dir, base_name)

    # unlink top and trim files if they already exist...
    # The calling function opens file in append mode "a", so important to start with empty file
    clean_up(top_file)
    
    # multicore stuff
    stdout_lock = Lock()
    top_file_lock = Lock()
    
    # Run octopus
    print("Execute Octopus on %s")
    run_multicore_function(seqbuddy.records, octopus, ["/usr/local/blastdbs/species_protein/Hydra_magnipapillata",
                                                       temporary_directory])

    print(temporary_directory.subdirs)
    # clean_up(["%s/PSSM_PRF_FILES" % out_dir, "%s/RAW_PRF_FILES" % out_dir, "%s/exec_times-bloctopus.txt" % out_dir,
    #          "%s/exec_times-octopus.txt" % out_dir])