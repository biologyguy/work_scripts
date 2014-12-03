#!/usr/bin/env python3
from MyFuncs import TempFile, TempDir, run_multicore_function
from Bio import SeqIO
from subprocess import Popen, PIPE
from os import remove
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
import datetime


class Prosite():
    def __init__(self, sequence):  # Sequence is a protein SeqRecord object
        self.tmp_dir = TempDir()
        self.job_id = ""
        sequence.seq.alphabet = IUPAC.IUPACProtein()
        self.sequence = sequence
        self.outfile = ""

    def run_prosite(self, client_path):
        tmp_file = TempFile()
        tmp_file.write(str(self.sequence.seq))
        file_location = "%s/%s" % (self.tmp_dir, str(tmp_file).split("/")[-1])
        tmp_file.save(file_location)

        output = Popen("%s --email biologyguy@gmail.com --outfile '%s/%s' --outputLevel 1 %s"
                       % (client_path, self.tmp_dir, self.sequence.id, file_location), shell=True,
                       stdout=PIPE).communicate()[0].decode()

        self.job_id = output.split("\n")[0]

        with open("%s/%s.out.txt" % (self.tmp_dir, self.sequence.id), "r") as in_file:
            self.outfile = in_file.read()

        for feature in self.outfile.split(">")[1:]:
            feat_type = re.match('EMBOSS_001 : (.*)', feature)
            feat_type = feat_type.groups(0)[0].split(" ")[1]
            feat_type = feat_type[:15]  # Need to limit the feature length, because gb format breaks otherwise
            spans = re.findall('([0-9]+ \- [0-9]+)', feature)
            for span in spans:
                span = span.split(" ")
                feature = SeqFeature(FeatureLocation(int(span[0]), int(span[2])), type=feat_type)
                self.sequence.features.append(feature)

        return self.sequence

    def save_prosite_file(self, out_path):
        today = datetime.date.today()
        today = "%s %s %s" % (today.strftime("%b"), today.day, today.year)
        with open(out_path, "w") as out_file:
            out_file.write("# Job ID: %s\t%s\n%s" % (self.job_id, today, self.outfile))


def run_interproscan(sequence, interpro_output_dir):  # This had not been fully implemented...
    tmp_file = TempFile()
    sequence.id = sequence.id
    tmp_file.write(str(sequence.seq))
    tmp_file.close()
    output = Popen("%s --email biologyguy@gmail.com --outfile '%s/%s' --outputLevel 1 --service interpro %s"
                   % (prosite_scan_client, interpro_output_dir, sequence.id, tmp_file.file), shell=True,
                   stdout=PIPE).communicate()[0].decode()

    project_id = output.split("\n")[0]

    with open("%s/%s.tsv.txt" % (interpro_output_dir, sequence.id), "r") as in_file:
        content = in_file.read()

    with open("%s/%s.tsv.txt" % (interpro_output_dir, sequence.id), "w") as out_file:
        out_file.write("# Job ID: %s\n%s" % (project_id, content))

    remove("%s/%s.gff.txt" % (interpro_output_dir, sequence.id))
    remove("%s/%s.htmltarball.html.tar.gz" % (interpro_output_dir, sequence.id))
    remove("%s/%s.log.txt" % (interpro_output_dir, sequence.id))
    remove("%s/%s.out.txt" % (interpro_output_dir, sequence.id))
    remove("%s/%s.sequence.txt" % (interpro_output_dir, sequence.id))
    remove("%s/%s.svg.svg" % (interpro_output_dir, sequence.id))
    remove("%s/%s.xml.xml" % (interpro_output_dir, sequence.id))


if __name__ == '__main__':
    import argparse
    import os
    from multiprocessing import Lock

    parser = argparse.ArgumentParser(prog="run_prosite.py", description="Converts a fasta file full of sequences into a"
                                                                        " genbank file full of sequences annotated with"
                                                                        " prosite motifs.")

    parser.add_argument("fasta_file", help="Where are the sequences you want analized", action="store")
    parser.add_argument("gb_file", help="Where would you like the output saved?", action="store")
    parser.add_argument("prosite_client", help="Location of ps_scan_py3.py", action="store")

    in_args = parser.parse_args()

    fasta_file = in_args.fasta_file  # "/Users/bondsr/Documents/Work/Innexin_evolution/Ctenophores/Ctenos_pep_fulls.fa"
    gb_file = os.path.abspath(in_args.gb_file)  # "/Users/bondsr/Documents/Work/Innexin_evolution/prosite_scan_files"
    prosite_scan_client = in_args.prosite_client  # "/Users/bondsr/Documents/public_scripts/ps_scan_py3.py"

    lock = Lock()
    with open(gb_file, "w") as ofile:
        ofile.truncate()

    def run_prosite(sequence):
        prosite = Prosite(sequence)
        sequence = prosite.run_prosite(prosite_scan_client)
        with lock:
            with open(gb_file, "a") as out_file:
                SeqIO.write(sequence, out_file, "gb")

    with open("%s" % fasta_file, "r") as ifile:
        sequences = list(SeqIO.parse(ifile, 'fasta'))

    run_multicore_function(sequences, run_prosite)