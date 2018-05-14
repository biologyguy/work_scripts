#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from buddysuite import SeqBuddy as Sb
from buddysuite import buddy_resources as br
import argparse
import re
from subprocess import Popen, PIPE
from multiprocessing import Lock, cpu_count
from Bio import SeqIO
from math import floor
from time import strftime
from shutil import move
from time import time


def mc_run_blast(records, args):
    blastdbs, evalue, threads = args

    tmp_file = br.TempFile()
    with open(tmp_file.path, "w") as _ofile:
        SeqIO.write(records, _ofile, "fasta")

    for blastdb in blastdbs:
        _cmd = "blastp -query %s -db %s -evalue %s -max_target_seqs 1000 -num_threads %s -dbsize 1000000000 " \
               "-outfmt '6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore " \
               "qlen slen nident'" % (tmp_file.path, blastdb, evalue, threads)

        _output = Popen(_cmd, stdout=PIPE, shell=True).communicate()
        _output = _output[0].decode()

        with file_lock:
            with open("%s/temp_blast_hits.csv" % in_args.outdir, "a") as _ofile:
                _ofile.write(_output)

    return


def backup():
    back_dir = "%s/%s_backup" % (in_args.outdir, strftime("%m_%d_%y_%H_%M_%s"))
    os.makedirs(back_dir)
    if os.path.isfile("%s/all_by_all.csv" % in_args.outdir):
        move("%s/all_by_all.csv" % in_args.outdir, "%s/all_by_all.csv" % back_dir)

    if os.path.isdir("%s/blastdbs" % in_args.outdir):
        move("%s/blastdbs" % in_args.outdir, "%s/blastdbs" % back_dir)


def prepare_previous():
    _root, _dirs, _files = next(os.walk("%s/blastdbs" % in_args.outdir))
    for db in _files:
        if db[-3:] == "phr":
            db = db.split("/")[-1]
            db = db[:-4]
            prev_blast_dbs.append("%s/blastdbs/%s" % (in_args.outdir, db))

    if os.path.isfile("%s/blastdbs/hash_map.csv" % in_args.outdir):
        with open("%s/blastdbs/hash_map.csv" % in_args.outdir, "r") as _ifile:
            for _line in _ifile:
                _next_hash, _name = _line.split(",")
                total_hash_map[_next_hash] = _name.strip()
                reverse_hash_map[_name.strip()] = _next_hash
    return


parser = argparse.ArgumentParser(prog='all-by-all_blast', description='Compare all sequences with blast2seq')

parser.add_argument('indir', help='Location of all files to be checked', action='store')
parser.add_argument('extensions', help="List each extension type to be checked", nargs="+", action='store')
parser.add_argument('--outdir', '-o', help='where should everything be saved.', default=os.getcwd())
parser.add_argument('--num_threads', '-nt', type=int, action='store', default=0,
                    help='Specify how many cores to dedicate to the job')
parser.add_argument('--e_value', '-E', action='store', type=float, default=0.00001,
                    help='Specify the E-value threshold')
parser.add_argument('--from_scratch', '-fs', action='store_true',
                    help="Force completely new search even if previous run is detected")
parser.add_argument('--original_names', '-on', action="store_true",
                    help='Do not change the names in input files to organism@seq_id. Just keep whatever is present.')
parser.add_argument('--strip_species', "-ss", help="Remove")
parser.add_argument('--log', '-l', help="Convert dynamic print to line-by-line logging", action="store_true")
parser.add_argument('--quiet', '-q', help="Suppress output to any DynamicPrint objects.", action="store_true")

in_args = parser.parse_args()

in_args.extensions = [re.sub("\.", "", ext) for ext in in_args.extensions]

print("*** Preparing to run all-by-all blastp ***")
# ### Prepare the script environment ### #
assert os.path.isdir(in_args.indir)
in_args.indir = os.path.abspath(in_args.indir)

assert os.path.isdir(in_args.outdir)
in_args.outdir = os.path.abspath(in_args.outdir)

cpus = floor((cpu_count() - 1))
cpus = 1 if cpus < 1 else cpus

if in_args.num_threads == 0:
    in_args.num_threads = cpus

if in_args.num_threads > cpus:
    in_args.num_threads = cpus

if os.path.isfile("%s/temp_blast_hits.csv" % in_args.outdir):
    os.remove("%s/temp_blast_hits.csv" % in_args.outdir)

# Don't delete anything if the -fs flag is thrown, just copy it all to a backup directory.
if in_args.from_scratch:
    backup()

root, dirs, files = next(os.walk(in_args.indir))
total_hash_map = {}
reverse_hash_map = {}
blast_dir = br.TempDir()
new_records_list = []
prev_records_list = []
prev_blast_dbs = []
new_blast_dbs = []
file_lock = Lock()
# ##################################### #


if in_args.strip_species:
    timer = br.RunTime(prefix="Run time: ")
    if os.path.isfile("%s/all_by_all.csv" % in_args.outdir):
        counter = 0
        if not in_args.quiet:
            print("Removing records for '%s' from all_by_all.csv" % in_args.strip_species)
            timer.start()
        infile = open("%s/all_by_all.csv" % in_args.outdir, "r")
        with open("%s/all_by_all_stripped.csv" % in_args.outdir, "w") as ofile:
            for line in infile:
                if line[0] == "#":
                    continue
                genes = line.split("\t")[:2]
                if not re.match("%s[^a-zA-Z0-9]" % in_args.strip_species, genes[0]) \
                        and not re.match("%s[^a-zA-Z0-9]" % in_args.strip_species, genes[1]):
                    ofile.write(line)
                else:
                    counter += 1
        infile.close()

        if not in_args.quiet:
            timer.end()
            print("Removed %s records\n" % counter)

    if os.path.isdir("%s/blastdbs" % in_args.outdir):
        if not in_args.quiet:
            print("Removing blast database\n")

        for root, dirs, files in br.walklevel("%s/blastdbs" % in_args.outdir):
            for _file in files:
                if re.match("%s\.p" % in_args.strip_species, _file):
                    os.remove("%s/%s" % (root, _file))

        if os.path.isfile("%s/blastdbs/hash_map.csv" % in_args.outdir):
            counter = 0
            if not in_args.quiet:
                print("Removing records for '%s' from hash_map.csv" % in_args.strip_species)
                timer.start()

            infile = open("%s/blastdbs/hash_map.csv" % in_args.outdir, "r")
            with open("%s/blastdbs/hash_map_stripped.csv" % in_args.outdir, "w") as ofile:
                for line in infile:
                    if line[0] == "#":
                        continue
                    gene = line.split(",")[1]
                    if not re.match("%s[^a-zA-Z0-9]" % in_args.strip_species, gene):
                        ofile.write(line)
                    else:
                        counter += 1

            infile.close()
            if not in_args.quiet:
                timer.end()
                print("Removed %s records\n" % counter)

    sys.exit()


if os.path.isdir("%s/blastdbs" % in_args.outdir) and os.path.isfile("%s/all_by_all.csv" % in_args.outdir):
    print("***Reading in previous data***")
    prepare_previous()
    print("    --> Done\n")
else:
    if not os.path.isdir("%s/blastdbs" % in_args.outdir):
        os.makedirs("%s/blastdbs" % in_args.outdir)

    with open("%s/all_by_all_headings.csv" % in_args.outdir, "w") as outfile:
        outfile.write("#query_id\tsubj_id\tperc_ident\talign_len\tmismatches\tgap_opens\tq_start\tq_end\ts_start\t"
                      "s_end\tevalue\tbit_score\tq_len\ts_len\n")

seq_files = []
for _file in files:
    extension = _file.split(".")[-1]
    name = re.sub(extension, "", _file.split("/")[-1])
    if extension in in_args.extensions and name not in prev_blast_dbs:
        seq_files.append("%s/%s" % (in_args.indir, _file))

print("***Hashing proteomes***")
for i in range(len(seq_files)):
    _file = seq_files[i]
    name = _file.split("/")[-1]
    name = "_".join(name.split(".")[:-1])

    seqbuddy = Sb.SeqBuddy(_file)
    seqbuddy = Sb.clean_seq(seqbuddy)

    if "%s/blastdbs/%s" % (in_args.outdir, name) in prev_blast_dbs:
        for record in seqbuddy.records:
            if not in_args.original_names:
                record.id = reverse_hash_map["%s@%s" % (name, record.id)]
            else:
                record.id = reverse_hash_map[record.id]
        prev_records_list += seqbuddy.records
        continue

    print(name)
    while True:
        redo = False
        seqbuddy = Sb.hash_ids(seqbuddy)
        for next_hash, line in seqbuddy.hash_map.items():
            # Re-hash if any duplicate hashes show up...
            hash_subset = []
            if next_hash in total_hash_map:
                for _hash in hash_subset:
                    del total_hash_map[_hash]
                print("Wow, duplicate hash detected. That's like winning the lottery! Let's make a new one...")
                redo = True
                break

            hash_subset.append(next_hash)
            if not in_args.original_names:
                total_hash_map[next_hash] = "%s@%s" % (name, line)
            else:
                total_hash_map[next_hash] = line
        if redo:
            continue

        seqbuddy = Sb.delete_metadata(seqbuddy)
        _file = "%s/%s" % (blast_dir.path, _file.split("/")[-1])
        new_records_list += seqbuddy.records

        with open(_file, "w") as ofile:
            SeqIO.write(seqbuddy.records, ofile, "fasta")

        seq_files[i] = _file
        break

with open("%s/blastdbs/hash_map.csv" % in_args.outdir, "w") as ofile:
    for _hash in total_hash_map:
        ofile.write("%s,%s\n" % (_hash, total_hash_map[_hash]))

print("\n***Creating new blast databases***")
for _file in seq_files:
    path = _file.split("/")
    name = ".".join(path[-1].split(".")[:-1])
    if "%s/blastdbs/%s" % (in_args.outdir, name) not in prev_blast_dbs:
        blast_db = "%s/blastdbs/%s" % (in_args.outdir, name)
        new_blast_dbs.append(blast_db)
        cmd = "makeblastdb -in %s -parse_seqids -dbtype prot -out %s -hash_index" % (_file, blast_db)
        Popen(cmd, shell=True).wait()

if not new_blast_dbs:
    sys.exit("No new proteomes detected. Exiting with nothing to be done.")

print("\n***Blasting %s new sequences against all databases***\n" % len(new_records_list))
if len(new_records_list) < in_args.num_threads * 1000:
    group_size = floor(len(new_records_list) / (in_args.num_threads - 1))
    group_size = 1 if group_size < 1 else group_size
    new_records_list = [new_records_list[i:i + group_size] for i in range(0, len(new_records_list), group_size)]
else:
    new_records_list = [new_records_list[i:i + 1000] for i in range(0, len(new_records_list), 1000)]

num_threads = 3 if len(new_records_list) <= cpus else 1
br.run_multicore_function(new_records_list, mc_run_blast, max_processes=in_args.num_threads, quiet=in_args.quiet,
                          func_args=[new_blast_dbs + prev_blast_dbs, in_args.e_value, num_threads], log=in_args.log)

if prev_records_list:
    print("\n***Blasting %s previous sequences against new databases***\n" % len(prev_records_list))
    if len(new_records_list) < in_args.num_threads * 1000:
        group_size = int(floor(len(prev_records_list) / (in_args.num_threads - 1)))
        prev_records_list = [prev_records_list[i:i + group_size] for i in range(0, len(prev_records_list), group_size)]
    else:
        prev_records_list = [prev_records_list[i:i + 1000] for i in range(0, len(prev_records_list), 1000)]

    num_threads = 3 if len(prev_records_list) <= in_args.num_threads else 1
    br.run_multicore_function(prev_records_list, mc_run_blast, [new_blast_dbs, in_args.e_value, num_threads],
                              max_processes=in_args.num_threads, quiet=in_args.quiet, log=in_args.log)

# Remove self-hits and convert any e-values of 0.0 to 1e-180
print("\n***Processing blast hits***")

blast_p_hits_handle = open("%s/temp_blast_hits.csv" % in_args.outdir, "r")

total_hits = 0
for total_hits, l in enumerate(blast_p_hits_handle):
    pass

total_hits += 1

# Convert seq IDs back from their hashes
blast_p_hits_handle.seek(0)
all_by_all_handle = open("%s/all_by_all.csv" % in_args.outdir, "a")

printer = br.DynamicPrint(quiet=in_args.quiet, log=in_args.log)
counter = 1
start_time = round(time())
for hit in blast_p_hits_handle:
    counter += 1
    if counter % 1000 == 0:
        printer.write("\t--> formatting %s of %s hits" % (counter, total_hits))
        printer.new_line()
    if hit == "":
        continue

    data = hit.split("\t")
    if data[0] == data[1]:
        continue

    if data[10] == '0.0':
        data[10] = '1e-181'

    data[0] = total_hash_map[data[0]]
    data[1] = total_hash_map[data[1]]
    all_by_all_handle.write("%s" % "\t".join(data))

printer.write("\t--> formatting %s of %s hits" % (counter, total_hits))

print("\n\tDone formatting in %s" % br.pretty_time(round(time()) - start_time))

os.remove("%s/temp_blast_hits.csv" % in_args.outdir)
