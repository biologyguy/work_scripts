#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shutil
import os
from rdmcl import helpers as hlp
import pandas as pd
from os.path import join
from subprocess import Popen
from buddysuite import buddy_resources as br
from buddysuite import SeqBuddy as Sb
import argparse

CPUS = br.usable_cpu_count()


def mc_psi_pred(seq_obj, args):
    outdir = args[0]
    if os.path.isfile(join(outdir, "%s.ss2" % seq_obj.id)):
        return
    result = run_psi_pred(seq_obj)
    with open(join(outdir, "%s.ss2" % seq_obj.id), "w") as ofile:
        ofile.write(result)
    return


def run_psi_pred(seq_rec):
    temp_dir = br.TempDir()
    pwd = os.getcwd()
    psipred_dir = join(hlp.SCRIPT_PATH, "psipred")
    os.chdir(temp_dir.path)
    with open("sequence.fa", "w") as ofile:
        ofile.write(seq_rec.format("fasta"))

    if shutil.which("psipred"):
        command = '''\
seq2mtx sequence.fa > {1}{3}{2}.mtx;
psipred {1}{3}{2}.mtx {0}{3}data{3}weights.dat {0}{3}data{3}weights.dat2 {0}{3}data{3}weights.dat3 > {1}{3}{2}.ss;
psipass2 {0}{3}data{3}weights_p2.dat 1 1.0 1.0 {1}{3}{2}.ss2 {1}{3}{2}.ss > {1}{3}{2}.horiz;
'''.format(psipred_dir, temp_dir.path, seq_rec.id, os.sep)

    else:
        data_weights = join(psipred_dir, "data", "weights")
        command = '''\
    {0}{3}bin{3}seq2mtx sequence.fa > {1}{3}{2}.mtx;
    {0}{3}bin{3}psipred {1}{3}{2}.mtx {4}.dat {4}.dat2 {4}.dat3 > {1}{3}{2}.ss;
    {0}{3}bin{3}psipass2 {4}_p2.dat 1 1.0 1.0 {1}{3}{2}.ss2 {1}{3}{2}.ss > {1}{3}{2}.horiz;
    '''.format(psipred_dir, temp_dir.path, seq_rec.id, os.sep, data_weights)

    Popen(command, shell=True).wait()
    os.chdir(pwd)
    with open(join(temp_dir.path, "%s.ss2" % seq_rec.id), "r") as ifile:
        result = ifile.read()
    return result


def read_ss2_file(path):
    ss_file = pd.read_csv(path, comment="#", header=None, delim_whitespace=True,
                          names=["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"])
    return ss_file


def main():
    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="homolog_tree_builder", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mRun PSI-PRED\033[m
  For Sofia, to do awesome stuff with

  Pass in a file of sequences, get secondary structure in return.
  
\033[1mUsage\033[m:
  run_psipred.py "/path/to/seqs" [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("seqs", help="Specify sequence file (most formats accepted)")
    positional.add_argument("save_ss2", action="store", help="Specify directory to save/read ss2 files.")
    
    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
    parser_flags.add_argument("-cpu", "--max_cpus", type=int, action="store", default=CPUS, metavar="",
                              help="Specify the maximum number of cores RD-MCL can use (default=%s)" % CPUS)

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version="1.0")
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()

    sequences = Sb.SeqBuddy(in_args.seqs)
    if not in_args.save_ss2:
        ss2_files = br.TempDir().path
    else:
        ss2_files = os.path.abspath(in_args.save_ss2)
        os.makedirs(ss2_files, exist_ok=True)

    br.run_multicore_function(sequences.records, mc_psi_pred, [ss2_files], max_processes=in_args.max_cpus)


if __name__ == '__main__':
    main()
