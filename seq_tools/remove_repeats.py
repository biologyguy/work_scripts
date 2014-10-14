#!/usr/bin/python3
from Bio import SeqIO, AlignIO, Seq
import argparse

parser = argparse.ArgumentParser(prog="Remove Repeats", description="Strip out replicate sequences from a fasta file and output to a new fasta file. Search on ID, sequence, or both",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', help='Location of input fasta file', action='store')
parser.add_argument('output', help='Location for output fasta file', action='store')
parser.add_argument('-s', '--search', help='Specify if you want to remove duplicate IDs, sequences, or both', choices=["id", "seq", "both"], action='store', default="id")

in_args = parser.parse_args()


def no_seq_repeat(seq_list, query, search_type):
    for seq in seq_list:
        if search_type in ["seq", "both"]:
            if str(query.seq) == str(seq.seq):
                return False
        if search_type in ["id", "both"]:
            if str(query.id) == str(seq.id):
                return False
    return True


with open(in_args.input, "r") as ifile:
    sequences = list(SeqIO.parse(ifile, "fasta"))

unique_seqs = []

for i in sequences:
    if no_seq_repeat(unique_seqs, i, in_args.search):
        unique_seqs.append(i)

with open(in_args.output, "w") as ofile:
    SeqIO.write(unique_seqs, ofile, "fasta")

print("Removed %s sequences from file, leaving %s sequences." % (len(sequences) - len(unique_seqs), len(unique_seqs)))
