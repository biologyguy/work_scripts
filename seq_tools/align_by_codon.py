#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC
import argparse
import sys
import re
from MyFuncs import *
from subprocess import Popen
from copy import copy

parser = argparse.ArgumentParser(prog="align_by_codon", description="Wrapper for mafft and muscle to do a translation"
                                                                    " alignment of CDSs.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('in_file', help="Path to fasta file. Only use in-frame CDSs, otherwise everything is going to be "
                                    "messed up.", action="store")
parser.add_argument('-p', '--params', help="Arguments that you want passed into command line (in quotes). "
                                           "Do not specify an out-file here, it will be removed if you try.",
                    action="store", default='')
parser.add_argument('-s', '--save_tmp', help="Keep the protein alignment file? Specify a directory.",
                    action="store", default=False)
parser.add_argument('-a', '--algorithm', help="Select the program you want to run", choices=["mafft", "muscle"],
                    action="store", default="mafft")

in_args = parser.parse_args()
in_file = os.path.abspath(in_args.in_file)

params = re.sub(">.*", "", in_args.params)


def guess_alphabet(fasta_file):
    valve = SafetyValve()
    with open(fasta_file, "r") as infile:
        seqs = SeqIO.parse(infile, "fasta")
        seq_concat = ""
        while len(seq_concat) < 1000:
            seq_concat += next(seqs).seq.upper()
            valve.test(seq_concat)

        percent_dna = float(seq_concat.count("A") + seq_concat.count("G") + seq_concat.count("T") +
                            seq_concat.count("C")) / float(len(seq_concat))
        if percent_dna > 0.9:
            return "nucl"
        else:
            return "prot"

# Do some sanity checking first
# Make sure we're working with DNA
if guess_alphabet(in_file) == "prot":
    sys.exit("You need to provide a DNA fasta file, not protein")

# Make sure the sequences are all divisible by three (most likely all codons)
with open(in_file, "r") as ifile:
    sequences = SeqIO.parse(ifile, "fasta")
    for seq in sequences:
        if len(seq.seq) % 3 != 0:
            sys.exit("%s in your fasta file is not a multiple of 3. Only provide full CDSs." % seq.id)


# Create protein translation
prot_file = TempFile()
prot_file.open()
with open(in_file, "r") as ifile:
    sequences = SeqIO.parse(ifile, "fasta")
    for seq in sequences:
        seq.seq = seq.seq.translate(to_stop=True)
        seq.alphabet = IUPAC.protein
        SeqIO.write(seq, prot_file.handle, "fasta")

prot_file.close()

out_file = TempFile()
print("Running alignment on translation", file=sys.stderr)

if in_args.algorithm == "mafft":
    Popen("mafft %s %s > %s" % (params, prot_file.path, out_file.path), shell=True).wait()

elif in_args.algorithm == "muscle":
    Popen("muscle -in %s -out %s %s" % (prot_file.path, out_file.path, params), shell=True).wait()

print("Converting translated alignment back to codons", file=sys.stderr)
out_file.open("r")
dna_file = open(in_file, "r")

if in_args.save_tmp:
    location = in_file.split("/")[-1].split(".")[:-1]
    out_file.save("%s/%s_aln_prot.fa" % (os.path.abspath(in_args.save_tmp), "_".join(location)))

alignment = next(AlignIO.parse(out_file.handle, "fasta",))
alignment = list(alignment)

sequences = SeqIO.to_dict(SeqIO.parse(dna_file, "fasta"))

new_alignment = []

for i in range(len(alignment)):
    # new_seq = Seq.Seq("", alphabet=IUPAC.unambiguous_dna)
    new_seq = copy(sequences[alignment[i].id])
    new_seq.seq = ""

    orig_dna = sequences[alignment[i].id].seq
    dna_pointer = 0

    align_prot = alignment[i].seq
    for aa in align_prot:
        if aa == "-":
            new_seq.seq += "---"
            continue
        else:
            new_seq.seq += orig_dna[dna_pointer:dna_pointer + 3]
            dna_pointer += 3

    new_alignment.append(new_seq)

out_file.close()
dna_file.close()

for seq in new_alignment:
    print(">%s\n%s\n" % (seq.id, seq.seq), file=sys.stdout)