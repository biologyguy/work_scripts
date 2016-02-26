#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014 

"""
Wraps a local install of the Octopus transmembrane prediction program (octopus.cbr.su.se).
Input is a sequence file (Biopython compliant), and each predicted TMD is appended to the seq feature list for printing
"""
from os import makedirs
from os.path import abspath
import re
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from MyFuncs import TempDir, run_multicore_function, walklevel
import SeqBuddy as Sb
from subprocess import Popen
from copy import copy


class AnnoTMD():
    """Takes a Bio.SeqRecord and runs octopus to find/annotate TMDs"""
    def __init__(self, sequence, blastdb, debug=False):
        sb = Sb.SeqBuddy([sequence])

        sequence.seq.alphabet = sb.alpha

        if sb.alpha != IUPAC.protein:
            dna_seq = copy(sequence)
            prot_seq = Seq(re.sub("\*", "", str(sequence.seq.translate())), alphabet=IUPAC.protein)
            sequence = SeqRecord(prot_seq, id=sequence.id, name=sequence.name, description=sequence.description)
        else:
            dna_seq = None

        sequence.seq = Seq(re.sub("\*", "", str(sequence.seq)), alphabet=sequence.seq.alphabet)

        tmp_dir = TempDir()
        out_dir = "%s/output" % tmp_dir.path
        makedirs(out_dir)

        with open("%s/%s.fa" % (tmp_dir.path, sequence.id), "w") as ofile:
            SeqIO.write(sequence, ofile, "fasta")

        with open("%s/NameFile.txt" % tmp_dir.path, "w") as file:
            file.write("%s\n" % sequence.id)

        popen_output = "> /dev/null 2>&1" if not debug else ""
        Popen("bloctopus %s/NameFile.txt %s %s /usr/local/bin/blastall /usr/local/bin/blastpgp %s "
              "/usr/local/bin/makemat -P%s" %
              (tmp_dir.path, tmp_dir.path, out_dir, blastdb, popen_output), shell=True).wait()

        Popen("octopus %s/NameFile.txt %s/PSSM_PRF_FILES %s/RAW_PRF_FILES %s -N%s"
              % (tmp_dir.path, out_dir, out_dir, out_dir, popen_output), shell=True).wait()

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
        out_file = abspath(out_file)
        with open(out_file, "w") as ofile:
            SeqIO.write(self.sequence, ofile, "gb")

    def __str__(self):
        _output = ""
        for i in range(len(self.TMDs)):
            _output += "TMD%s\t%s-%s\n" % (i + 1, self.TMDs[i][0], self.TMDs[i][1])
        return _output.strip()


if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog="annotate_TMDs.py", description="Identify transmembrane domains")

    parser.add_argument("in_file", help="Location of DNA or protein input file", action="store")
    parser.add_argument("blastdb", help="Path to BLASTP database used by Octopus", action="store")
    parser.add_argument("-i", "--in_place", help="Overwrite the original file. BE CAREFUL!", action="store_true")

    in_args = parser.parse_args()

    def get_tmds(seq_rec, args):
        _temp_dir = args[0]
        tmds = AnnoTMD(seq_rec, in_args.blastdb)
        tmds.save("%s/%s.gb" % (_temp_dir, seq_rec.id))

    seqs = Sb.SeqBuddy(in_args.in_file)
    temp_dir = TempDir()
    run_multicore_function(seqs.records, get_tmds, [temp_dir.path], out_type=sys.stderr)

    output = ""

    root, dirs, files = next(walklevel(temp_dir.path))

    for file in files:
        with open("%s/%s" % (root, file), "r") as ifile:
            output += "%s\n" % ifile.read()

    if in_args.in_place:
        with open(in_args.in_file, "w") as ofile:
            ofile.write(output)
        print("File over-written at:\n%s" % in_args.in_file)
    else:
        print(output)
