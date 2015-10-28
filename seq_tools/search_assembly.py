#!/usr/bin/python3
# -*- coding: utf-8 -*-

from MyFuncs import *
import sys
import os
import shutil
from subprocess import Popen, PIPE
from multiprocessing import Lock
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from SeqBuddy import blast, SeqBuddy, _guess_alphabet, rename


class SeqTools():
    def __init__(self, fasta_file, outdir, seq_type=False):
        self.fasta_file = os.path.abspath(fasta_file)

        if _guess_alphabet(SeqBuddy(self.fasta_file)) == "prot":
            sys.exit("Error: Can't run transdecoder on a protein database.")

        self.name = "_".join(fasta_file.split("/")[-1].split(".")[:-1])
        self.outdir = os.path.abspath(outdir)

        if not os.path.exists(self.outdir):
            sys.exit("Error: The indicated output directory does not exist. Please create it first.")

        self.temp_dir = TempDir()
        self.trans_d_dir = "%s/transdecoder/%s" % (self.outdir, self.name)
        self.blastdb_dir = "%s/blastdbs" % self.outdir
        self.blast_out_dir = "%s/blast_out" % self.outdir

        self.transd_run = False
        self.makeblastdb_run = False
        self.blast_run = False

        self.lock = Lock()

    def transdecoder(self, save_tmp=False, _force=False):
        self.transd_run = True
        if os.path.exists(self.trans_d_dir):
            check_extensions = []
            pwd, dirs, files = next(walklevel(self.trans_d_dir))
            for next_file in files:
                check_extensions.append(next_file.split(".")[-1])

            if check_extensions == ["bed", "cds", "gff3", "mRNA", "pep"] and not _force:
                print("Transdecoder output detected for %s, using previous data. (You can override by force)" % self.name)
                return False

        cwd = os.getcwd()
        os.chdir(self.temp_dir.path)

        if not os.path.exists("%s/transdecoder" % self.outdir):
            print("Creating %s/transdecoder" % self.outdir)
            os.mkdir("%s/transdecoder" % self.outdir)

        if not os.path.exists("%s" % self.trans_d_dir):
            print("Creating %s" % self.trans_d_dir)
            os.mkdir(self.trans_d_dir)

        print("Running TransDecoder on %s\ntransdecoder --CPU 21 -t %s\n" % (self.name, self.fasta_file))
        Popen("transdecoder --CPU 21 -t %s" % self.fasta_file, shell=True).wait()

        pwd, dirs, files = next(walklevel("./"))
        print(dirs)
        for next_file in files:
            if next_file.split(".")[0] == self.name:
                extension = next_file.split(".")[-1]
                new_name = "%s_transD.%s" % (self.name, extension)
                shutil.move(next_file, "%s/%s" % (self.trans_d_dir, new_name))
                # Clean up the ids of the output files
                if extension in ["cds", "pep", "mRNA"]:
                    out_seqs = SeqBuddy("%s/%s" % (self.trans_d_dir, new_name))
                    rename(out_seqs, ".*\|m\.", "m")
                    with open("%s/%s" % (self.trans_d_dir, new_name), "w") as ofile:
                        SeqIO.write(out_seqs, ofile, "fasta")

        if save_tmp:
            temp_files = dirs[0]
            pwd, dirs, files = next(walklevel(temp_files))
            if not os.path.exists("%s/temp_files" % self.trans_d_dir):
                os.mkdir("%s/temp_files" % self.trans_d_dir)

            for next_file in files:
                shutil.move("%s/%s" % (temp_files, next_file), "%s/temp_files/%s" % (self.trans_d_dir, next_file))

        print("\nTransDecoder done, output files in %s" % self.trans_d_dir)
        os.chdir(cwd)
        return True

    def makeblastdb(self, _force=False):
        self.makeblastdb_run = True

        if not os.path.exists(self.blastdb_dir):
            os.mkdir(self.blastdb_dir)

        # Make both prot and nucl databases
        if not os.path.exists("%s/%s_prot.phr" % (self.blastdb_dir, self.name)) or _force:
            print("Running makeblastdb on %s" % self.name)
            print("makeblastdb -in %s/%s_transD.pep -dbtype prot -out %s/%s_prot -parse_seqids\n" %
                  (self.trans_d_dir, self.name, self.blastdb_dir, self.name))
            Popen("makeblastdb -in %s/%s_transD.pep -dbtype prot -out %s/%s_prot -parse_seqids" %
                  (self.trans_d_dir, self.name, self.blastdb_dir, self.name), shell=True).wait()

        else:
            print("Protein blast db detected for %s, using previous data. (You can override with -f)" % self.name)
            return False

        if not os.path.exists("%s/%s_prot.nhr" % (self.blastdb_dir, self.name)) or _force:
            print("Running makeblastdb on %s" % self.name)
            print("makeblastdb -in %s/%s_transD.cds -dbtype nucl -out %s/%s_nucl -parse_seqids\n" %
                  (self.trans_d_dir, self.name, self.blastdb_dir, self.name))
            Popen("makeblastdb -in %s/%s_transD.cds -dbtype nucl -out %s/%s_nucl -parse_seqids" %
                  (self.trans_d_dir, self.name, self.blastdb_dir, self.name), shell=True).wait()

        else:
            print("CDS blast db detected for %s, using previous data. (You can override with -f)" % self.name)
            return False

        return True

    def blast(self, query_file, _force=False, evalue=0.01):
        if not os.path.exists(query_file):
            sys.exit("Error: The indicated query file does not exist. I need some sequences to search with.")

        query_file = os.path.abspath(query_file)
        query_name = query_file.split("/")[-1].split(".")[0]

        if not self.makeblastdb_run:
            self.makeblastdb()

        if not os.path.exists(self.blast_out_dir):
            os.mkdir(self.blast_out_dir)

        query_seqs = SeqBuddy(query_file)
        seq_type = _guess_alphabet(query_seqs)

        blast_hits = blast(query_seqs, "%s/%s_%s" % (self.blastdb_dir, self.name, seq_type))

        with open("%s/%s|%s_query.fa" % (self.blast_out_dir, self.name, query_name), "w") as ofile:
            for seq in blast_hits.seqs:
                SeqIO.write(seq, ofile, "fasta")

        self.blast_run = True
        return True

    def _parse_hits(self, hit_id, args):
        query_name = args[0]
        seq = Popen("pull_fasta_rec %s/%s_transD.cds '%s'" % (self.trans_d_dir, self.name, hit_id),
                    shell=True, stdout=PIPE).communicate()[0]

        temp_file = TempFile()
        with open(temp_file.path, "w") as ofile:
            ofile.write(seq.decode("utf-8"))

        with open(temp_file.path, "r") as ifile:
            seq = SeqIO.read(ifile, 'fasta')

        with self.lock:
            with open("%s/%s|%s_query.fa" % (self.blast_out_dir, self.name, query_name), "a") as ofile:
                SeqIO.write(seq, ofile, "fasta")

    def _get_records(self, query_name):
        with open("%s/%s|%s_query.txt" % (self.blast_out_dir, self.name, query_name), "r") as ifile:
            blast_results = SearchIO.parse(ifile, "blast-tab")
            records = list(blast_results)
        return records

    def _guess_alphabet(self):
        valve = SafetyValve()
        with open(self.fasta_file, "r") as ifile:
            sequences = SeqIO.parse(ifile, "fasta")
            seq_concat = ""
            while len(seq_concat) < 1000:
                seq_concat += next(sequences).seq.upper()
                valve.test(seq_concat)

            percent_dna = float(seq_concat.count("A") + seq_concat.count("G") + seq_concat.count("T") + seq_concat.count("C")) / float(len(seq_concat))
            if percent_dna > 0.9:
                return "nucl"
            else:
                return "prot"


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="sequence_tools",
                                     description="Wrapper tools for Transdecoder and some Blast++ stuff")

    parser.add_argument('input_file', help='FASTA file of DNA or Protein sequences')

    parser.add_argument('-o', '--out_dir',
                        help='Specify the directory where you want your output files. Default is current working dir.',
                        action='store', default=os.getcwd())
    parser.add_argument('-t', '--seq_type',
                        help='Specify whether sequences are protein or dna, otherwise the program will guess.',
                        type=str, choices=["prot", "nucl"], default=False)
    parser.add_argument('-d', '--transdecoder', help='Run transdecoder on input sequences.', action='store_true',
                        default=False)
    parser.add_argument('-m', '--makeblastdb',
                        help='Run makeblastdb on input sequences, or transdecoder file if that was run first.',
                        action='store_true', default=False)
    parser.add_argument('-b', '--blast',
                        help="Run blastp or blastn, and needs an input query file (defaults to blastp if transdecoder "
                             "has been run, but can be over-ridden by specifying 'prot' or 'nucl' with --seq_type)",
                        action='store')
    parser.add_argument('-f', '--force', nargs='+', choices=["transD", "blast", "makeblastdb"],
                        help="Indicate if you want previous runs to be over-written", default=[])

    in_args = parser.parse_args()

    transome = SeqTools(in_args.input_file, os.path.abspath(in_args.out_dir), seq_type=in_args.seq_type)

    print(transome.name)

    if in_args.transdecoder:
        force = True if ("transD" in in_args.force) else False
        transome.transdecoder(_force=force)

    if in_args.makeblastdb:
        force = True if ("makeblastdb" in in_args.force) else False
        transome.makeblastdb(_force=force)

    if in_args.blast:
        force = True if ("blast" in in_args.force) else False
        transome.blast(os.path.abspath(in_args.blast), _force=force)

