#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Nov 18 2014 

"""
DESCRIPTION OF PROGRAM
AlignmentBuddy is a general wrapper for popular DNA and protein alignment programs that handles format conversion
and allows maintencance of rich feature annotation following alignment.
"""

import argparse
import MyFuncs
from Bio import SeqIO
import os
import sys
import shutil


class NewClass():
    """DESCRIPTION OF CLASS"""
    def __init__(self):
        self.x = 1

    def class_def(self):
        self.x = 1
        return self.x


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


def def2():
    """DESCRIPTION OF FUNC"""
    x = 1
    return x


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="alignBuddy", description="Sequience alignment with a splash of Kava")

    parser.add_argument("sequence_file", help="Where are the sequences you want to align?", action="store")
    parser.add_argument("alignment_package", help="Pick your aligner. If the package is not in your $PATH, "
                                                  "specify the location with -a",
                        choices=["mafft", "prank", "pagan", "muscle", "clustalw"])
    parser.add_argument('-p', '--parameters', help="Arguments that you want passed to the alignment software (in quotes). "
                                                   "Do not specify an out-file here, it will be removed if you try.",
                        action="store", default='')
    parser.add_argument("-a", "--align_binary", help="Specify the path to your alignment package if not in $PATH",
                        action="store")
    parser.add_argument("-if", "--in_format", help="What format is your sequence in?",
                        choices=["fasta", "stockholm", "phylip", "nexus", "gb", "genbank"], action="store",
                        type=str, default='genbank')
    parser.add_argument("-of", "--out_format", help="How would you like the results returned?",
                        choices=["fasta", "stockholm", "phylip", "phylip-relaxed", "nexus", "clustal"], action="store",
                        type=str, default='stockholm')

    #parser.add_argument("-c", "--choice", help="", type=str, choices=["", ""], default=False)
    #parser.add_argument("-m", "--multi_arg", nargs="+", help="", default=[])

    in_args = parser.parse_args()
    if not in_args.align_binary:
        if not shutil.which(in_args.alignment_package):
            sys.exit("Error: Unable to locate %s in your $PATH. Please specify a path to the binary with the -a flag."
                     % in_args.alignment_package)

    else:
        align_binary = os.path.abspath(in_args.align_binary)
        if not os.path.isfile(align_binary):
            sys.exit("Error: Unable to resolve the provided path to %s.\nUser input: %s" %
                     (in_args.alignment_package, align_binary))

    seq_file = os.path.abspath(in_args.sequence_file)
    if not os.path.isfile(seq_file):
        sys.exit("Error: Unable to find the sequence file specified. Is the path correct?")

    with open(seq_file, "r") as ifile:
        in_seqs = SeqIO.to_dict(SeqIO.parse(ifile, in_args.in_format))
        print(in_seqs["AAEL006726"])
