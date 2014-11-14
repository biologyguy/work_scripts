#!/usr/bin/python3
from MyFuncs import TempFile, run_multicore_function
from Bio import SeqIO
from subprocess import Popen, PIPE
from os import remove
import argparse

output_dir = "/Users/bondsr/Documents/Work/Innexin_evolution/prosite_scan_files"
infile = "/Users/bondsr/Documents/Work/Innexin_evolution/all_innexins_pep.fa"
prosite_scan_client = "/Users/bondsr/Documents/public_scripts/ps_scan_py3.py"


def run_prosite(sequence):
    tmp_file = TempFile()
    sequence.id = sequence.id.split("|")[0]
    with open(tmp_file.file, "w") as ofile:
        ofile.write(str(sequence.seq))

    output = Popen("%s --email biologyguy@gmail.com --outfile '%s/%s' --outputLevel 1 %s"
                   % (prosite_scan_client, output_dir, sequence.id, tmp_file.file), shell=True, stdout=PIPE).communicate()[0].decode()

    project_id = output.split("\n")[0]

    with open("%s/%s.out.txt" % (output_dir, sequence.id), "r") as ifile:
        content = ifile.read()

    with open("%s/%s.out.txt" % (output_dir, sequence.id), "w") as ofile:
        ofile.write("# Job ID: %s\n%s" % (project_id, content))

    remove("%s/%s.sequence.txt" % (output_dir, sequence.id))


with open("%s" % infile, "r") as ifile:
    sequences = list(SeqIO.parse(ifile, 'fasta'))

run_multicore_function(sequences, run_prosite)
