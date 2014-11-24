#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 

"""
DESCRIPTION OF PROGRAM
"""
import os
import re
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from MyFuncs import TempDir, run_multicore_function
import seq_tools
from subprocess import Popen
from copy import copy


class AnnoTMD():
    """Takes a Bio.SeqRecord and runs octopus to find/annotate TMDs"""
    def __init__(self, sequence, blastdb, debug=False):
        alpha_guess = seq_tools.guess_alphabet(str(sequence.seq))
        sequence.seq.alphabet = IUPAC.unambiguous_dna if alpha_guess == "nucl" \
            else IUPAC.protein

        if alpha_guess == "nucl":
            dna_seq = copy(sequence)
            prot_seq = Seq(re.sub("\*", "", sequence.seq.translate()), alphabet=IUPAC.protein)
            sequence = SeqRecord(prot_seq, id=sequence.id)
        else:
            dna_seq = None

        sequence.seq = Seq(re.sub("\*", "", str(sequence.seq)), alphabet=sequence.seq.alphabet)

        tmp_dir = TempDir()
        out_dir = "%s/output" % tmp_dir.dir
        os.makedirs(out_dir)

        with open("%s/%s.fa" % (tmp_dir.dir, sequence.id), "w") as ofile:
            SeqIO.write(sequence, ofile, "fasta")

        with open("%s/NameFile.txt" % tmp_dir.dir, "w") as file:
            file.write("%s\n" % sequence.id)

        popen_output = "> /dev/null 2>&1" if not debug else ""
        Popen("bloctopus %s/NameFile.txt %s %s /usr/local/bin/blastall /usr/local/bin/blastpgp %s "
              "/usr/local/bin/makemat -P%s" %
              (tmp_dir.dir, tmp_dir.dir, out_dir, blastdb, popen_output), shell=True).wait()

        Popen("octopus %s/NameFile.txt %s/PSSM_PRF_FILES %s/RAW_PRF_FILES %s -N%s"
              % (tmp_dir.dir, out_dir, out_dir, out_dir, popen_output), shell=True).wait()

        with open("%s/%s.top" % (out_dir, sequence.id), "r") as ifile:
            top_file = SeqIO.read(ifile, "fasta", alphabet=IUPAC.protein)

        tmds = []
        counter = 1
        for membrane in re.finditer("M+", str(top_file.seq)):
            tmds.append(membrane.span())
            feature = SeqFeature(location=FeatureLocation(tmds[-1][0], tmds[-1][1]), type="TMD%s" % counter)
            top_file.features.append(feature)
            sequence.features.append(feature)
            counter += 1

        # Accessible variables
        self.sequence = sequence
        self.dna_sequence = dna_seq
        self.TMDs = tmds
        self.top_file = top_file

    def save(self, out_file):
        out_file = os.path.abspath(out_file)
        with open(out_file, "w") as ofile:
            SeqIO.write(self.sequence, ofile, "gb")

    def __str__(self):
        output = ""
        for i in range(len(self.TMDs)):
            output += "TMD%s\t%s-%s\n" % (i + 1, self.TMDs[i][0], self.TMDs[i][1])
        return output.strip()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="annotate_TMDs.py", description="Identify all the transmembrane domains in "
                                                                          "a fasta file, and send the output somewhere")

    parser.add_argument("in_file", help="Location of DNA or protein fasta file", action="store")
    parser.add_argument("out_dir", help="Where do you want the new genbank files sent", action="store")
    parser.add_argument("blastdb", help="Path to BLASTP database used by Octopus", action="store")

    in_args = parser.parse_args()

    def get_tmds(seq_rec):
        tmds = AnnoTMD(seq_rec, in_args.blastdb)
        tmds.save("%s/%s.gb" % (in_args.out_dir, seq_rec.id))

    with open(in_args.in_file, "r") as ifile:
        seqs = list(SeqIO.parse(ifile, "fasta"))

    run_multicore_function(seqs, get_tmds)