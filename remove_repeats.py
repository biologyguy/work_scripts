#!/usr/bin/python
from Bio import SeqIO, AlignIO, Seq
import argparse

parser = argparse.ArgumentParser(prog="Remove Repeats", description="Strip out replicate sequences from a fasta file and output to a new fasta file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', help='Location of input fasta file', action='store')
parser.add_argument('output', help='Location for output fasta file', action='store')

in_args = parser.parse_args()


def no_repeat(seq_list, query):
    for i in seq_list:
        if str(query) == str(i.seq):
            return False
    return True


with open(in_args.input, "r") as ifile:
    sequences = list(SeqIO.parse(ifile, "fasta"))

unique_seqs = []
for i in sequences:
    if no_repeat(unique_seqs, i.seq):
        unique_seqs.append(i)

print("Removed %s sequences from file, leaving %s sequences." % (len(sequences)-len(unique_seqs), len(unique_seqs)))

with open(in_args.output, "w") as ofile:
    SeqIO.write(unique_seqs, ofile, "fasta")