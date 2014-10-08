#!/usr/bin/python

from Bio import SeqIO

with open("Maur_transcriptome_v1.fasta", "r") as ifile:
    seqs = SeqIO.parse(ifile, "fasta")

    with open("Maur_transcritpome_v1_short.fasta", "w") as ofile:
        for seq in seqs:
            seq.description = ""
            SeqIO.write(seq, ofile, "fasta")