#!/usr/bin/python3
# -*- coding: utf-8 -*-

from MyFuncs import *
import sys
import os
from subprocess import Popen, PIPE
from multiprocessing import Lock
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq


class SeqTools():
    def __init__(self, fasta_file, outdir, seq_type=False):
        self.fasta_file = os.path.abspath(fasta_file)
        self.name = "_".join(fasta_file.split("/")[-1].split(".")[:-1])
        self.outdir = os.path.abspath(outdir)

        if not os.path.exists(self.outdir):
            sys.exit("Error: The indicated output directory does not exist. Please create it first.")

        self.temp_dir = TempDir()
        self.trans_d_dir = "%s/transdecoder/%s" % (self.outdir, self.name)
        self.blastdb_dir = "%s/blastdbs" % self.outdir
        self.blast_out_dir = "%s/blast_out" % self.outdir

        self.seq_type = seq_type if (seq_type and seq_type in ["prot", "nucl"]) else self._guess_alphabet()

        self.transd_run = False
        self.makeblastdb_run = False
        self.blast_run = False

        self.lock = Lock()

    def transdecoder(self, save_tmp=False, force=False):
        if self.seq_type == "prot":
            print("Can't run transdecoder on a protein database.")
            return False

        self.transd_run = True
        if os.path.exists(self.trans_d_dir):
            check_extensions = []
            pwd, dirs, files = next(walklevel(self.trans_d_dir))
            for next_file in files:
                check_extensions.append(next_file.split(".")[-1])

            if check_extensions == ["bed", "cds", "gff3", "mRNA", "pep"] and not force:
                print("Transdecoder output detected for %s, using previous data. (You can override by force)" % self.name)
                return False

        cwd = os.getcwd()
        os.chdir(self.temp_dir.dir)

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
                new_name = "%s_transD.%s" % (self.name, next_file.split(".")[-1])
                os.rename(next_file, "%s/%s" % (self.trans_d_dir, new_name))

        if save_tmp:
            temp_files = dirs[0]
            pwd, dirs, files = next(walklevel(temp_files))
            if not os.path.exists("%s/temp_files" % self.trans_d_dir):
                os.mkdir("%s/temp_files" % self.trans_d_dir)

            for next_file in files:
                os.rename("%s/%s" % (temp_files, next_file), "%s/temp_files/%s" % (self.trans_d_dir, next_file))

        print("\nTransDecoder done, output files in %s" % self.trans_d_dir)
        os.chdir(cwd)
        return True

    def makeblastdb(self, force=False):
        self.makeblastdb_run = True

        if not os.path.exists(self.blastdb_dir):
            os.mkdir(self.blastdb_dir)

        if not os.path.exists("%s/%s_prot.phr" % (self.blastdb_dir, self.name)) or force:
            print("Running makeblastdb on %s" % self.name)
            print("makeblastdb -in %s -dbtype %s -out %s/%s_%s\n" % (self.fasta_file, self.seq_type, self.blastdb_dir, self.name, self.seq_type))
            Popen("makeblastdb -in %s -dbtype %s -out %s/%s_%s" % (self.fasta_file, self.seq_type, self.blastdb_dir, self.name, self.seq_type), shell=True).wait()

        else:
            print("Blast db detected for %s, using previous data. (You can override by force)" % self.name)
            return False

        return True

    def blast(self, query_file, force=False, keep_blast=False, evalue=0.0001):
        if not os.path.exists(query_file):
            sys.exit("Error: The indicated query file does not exist. I need some sequences to search with.")

        query_file = os.path.abspath(query_file)
        query_name = query_file.split("/")[-1].split(".")[0]

        if not self.makeblastdb_run:
            self.makeblastdb()

        if not os.path.exists(self.blast_out_dir):
            os.mkdir(self.blast_out_dir)

        blast_program = "blastp" if self.seq_type == "prot" else "blastn"
        if not os.path.exists("%s/%s|%s_query.txt" % (self.blast_out_dir, self.name, query_name)) or force:
            print("Running %s with %s on %s" % (blast_program, query_name, self.name))
            print("%s -db %s/%s_%s -query %s -out '%s/%s|%s_query.txt' -num_threads 20 -evalue %s -outfmt 6\n" % (blast_program, self.blastdb_dir, self.name, self.seq_type, query_file, self.blast_out_dir, self.name, query_name, evalue))
            Popen("%s -db %s/%s_%s -query %s -out '%s/%s|%s_query.txt' -num_threads 20 -evalue %s -outfmt 6" % (blast_program, self.blastdb_dir, self.name, self.seq_type, query_file, self.blast_out_dir, self.name, query_name, evalue), shell=True).wait()

        else:
            print("%s output file detected for %s and %s, using previous data. (You can override by force)\n" % (blast_program, self.name, query_name))
            return False

        print("Collecting blast hits into FASTA format\n%s\n" % "%s/%s|%s_query.fa" % (self.blast_out_dir, self.name, query_name))
        records = self._get_records(query_name)

        hit_ids = []
        for record in records:
            for hsp in record.hsps:
                hit_id = hsp.hit_id

                if hit_id in hit_ids:
                    continue

                hit_ids.append(hit_id)

        print("Collecting match records:")
        with open("%s/%s|%s_query.fa" % (self.blast_out_dir, self.name, query_name), "w") as ofile:
                ofile.write("")

        run_multicore_function(hit_ids, self._parse_hits, [query_name])

        self.blast_run = True
        return True

    def _parse_hits(self, hit_id, args):
        query_name = args[0]
        seq = Popen("pull_fasta_rec %s '%s'" % (self.fasta_file, hit_id),
                    shell=True, stdout=PIPE).communicate()[0]

        temp_file = TempFile()
        with open(temp_file.file, "w") as ofile:
            ofile.write(seq.decode("utf-8"))

        with open(temp_file.file, "r") as ifile:
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

    parser = argparse.ArgumentParser(prog="sequence_tools", description="Wrapper tools for Transdecoder and some Blast++ stuff")

    parser.add_argument('input_file', help='FASTA file of DNA or Protein sequences')

    parser.add_argument('-o', '--out_dir', help='Specify the directory where you want your output files. Default is current working dir.', action='store', default=os.getcwd())
    parser.add_argument('-t', '--seq_type', help='Specify whether sequences are protein or dna, otherwise the program will guess.', type=str, choices=["prot", "nucl"], default=False)
    parser.add_argument('-d', '--transdecoder', help='Run transdecoder on input sequences.', action='store_true', default=False)
    parser.add_argument('-m', '--makeblastdb', help='Run makeblastdb on input sequences, or transdecoder file if that was run first.', action='store_true', default=False)
    parser.add_argument('-b', '--blast', help="Run blastp or blastn, and needs an input query file (defaults to blastp if transdecoder has been run, but can be over-ridden by specifying 'prot' or 'nucl' with --seq_type)", action='store')
    parser.add_argument('-f', '--force', nargs='+', choices=["transD", "blast", "makeblastdb"], help="Indicate if you want previous runs to be over-written", default=[])

    in_args = parser.parse_args()

    transome = SeqTools(in_args.input_file, os.path.abspath(in_args.out_dir), seq_type=in_args.seq_type)

    print(transome.name)

    if in_args.transdecoder:
        force = True if ("transD" in in_args.force) else False
        transome.transdecoder(force=force)

    if in_args.makeblastdb:
        force = True if ("makeblastdb" in in_args.force) else False
        transome.makeblastdb(force=force)

    if in_args.blast:
        force = True if ("blast" in in_args.force) else False
        transome.blast(os.path.abspath(in_args.blast), force=force)

