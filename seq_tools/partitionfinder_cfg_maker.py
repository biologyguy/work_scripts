#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Oct 3 2014

"""
Creates .cgf configuration files for partitionfinder(protein), given an input file that specifies where alignment files can be
found and how to break them up into data blocks. The actual program can be run directly from here as well.
For some PartitionFinder references, see:   - Lanfear et al., 2012 doi: 10.1093/molbev/mss020
                                            - Lanfear et al., 2014 doi: 10.1186/1471-2148-14-82
"""

from Bio import AlignIO
import argparse
import os
import sys
import re
from subprocess import Popen
from time import clock
import MyFuncs
import seq_tools


def make_cfg(location, _blocks, phylip_file):
    _output = "## ALIGNMENT FILE ##\n"
    _output += "alignment = %s_hashed.phy;\n\n" % phylip_file

    _output += "## BRANCHLENGTHS: linked | unlinked ##\n"
    _output += "branchlengths = linked;\n\n"  # options linked/unlinked

    _output += "## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | beast | <list> ##\n"
    _output += "##              for PartitionFinderProtein: all_protein | <list> ##\n"
    models = "all" if in_args.type == "nucl" else "all_protein"
    _output += "models = %s;\n\n" % models  # options all/some group pf models to test

    _output += "# MODEL SELECCTION: AIC | AICc | BIC #\n"
    _output += "model_selection = AIC;\n\n"

    _output += "## DATA BLOCKS: see manual for how to define ##\n"
    _output += "[data_blocks]\n"

    for _block in _blocks:
        _block = _block.split(",")
        if in_args.codon_pos:
            for j in range(3):
                _output += "%s_pos%s = %s-%s\\3;\n" % (_block[0], j + 1, int(_block[1]) + j, _block[2])
        else:
            _output += "%s = %s-%s;\n" % (_block[0], _block[1], _block[2])

    _output += "\n## SCHEMES, search: all | greedy | rcluster | hcluster | user ##\n"
    _output += "[schemes]\n"
    _output += "search = %s;\n\n" % in_args.algorithm

    with open("%s/partition_finder.cfg" % location, "w") as cfg:
        cfg.write(_output)

    return

parser = argparse.ArgumentParser(prog="partitionfinder_cgf_maker",
                                 description="Creates a new .cgf configuration file for partitionfinder",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", '--blocks_file',
                    help='Location of data_blocks.csv file that lists the alignments and how to break them up.',
                    action='store')
parser.add_argument("-ff", "--file_format", help="Print an example 'data_blocks.csv' file and exit.",
                    action="store_true", default=False)
parser.add_argument("-o", "--out_dir", help='Where would you like the output directories and files?',
                    action='store', default=os.getcwd())
parser.add_argument("-a", "--algorithm",
                    help="Specify which PartitionFinder algorithm you want to run. Greedy is most accurate, but to "
                         "slow for more than 100 partitions, strict hierarchical clustering (hcluster) is fastest, but "
                         "pretty terrible, and relaxed heirachiacal clustering (rcluster) is a good mix of speed and "
                         "accuracy for > 100 partitions",
                    choices=["greedy", "rcluster", "hcluster"], default="greedy")
parser.add_argument("-c", "--codon_pos", help="Break every data block up into the three different codon positions",
                    action="store_true", default=False)
parser.add_argument("-t", "--type", help="DNA or protein?", action="store", choices=["prot", "nucl", "guess"],
                    default="guess")
parser.add_argument("-r", "--run_partfinder", help="If you want to call PartitionFinder right away, go for it!",
                    action="store_true", default=False)
parser.add_argument("-f", "--force", help="If you really want to force PartitionFinder to run a job that already ran.",
                    action="store_true", default=False)

in_args = parser.parse_args()

start_time = round(clock())

if in_args.file_format:
    output = "# This is a comment and will be ignored\n"
    output += "Cnidaria_einsi.fa\t# alignment file name\n"
    output += "N-term,1,171\t# block name and start-stop positions. Names must be unique to this set of blocks, and " \
              "spaces will be replaced by underscores\n"
    output += "TMD1,172,237\n"
    output += "ECL1,238,585\n"
    output += "TMD2,586,651\n"
    output += "ICL,652,924\n"
    output += "TMD3,925,993\n"
    output += "ECL2,994,1278\n"
    output += "TMD4,1279,1350\n"
    output += "C-term,1351,1866\n"
    output += "//\t\t# End the set if more than one file is to be analyzed.\n"
    output += "Cnidaria_einsi_trimal.fasta\t# The file names are going to be used to make new directories\n"
    output += "N-term,1,75\n"
    output += "TMD2,76,141\n"
    output += "ECL2,142,390\n"
    output += "TMD3,391,456\n"
    output += "ICL,457,588\n"
    output += "TMD4,589,657\n"
    output += "ECL3,658,894\n"
    output += "TMD5,895,966\n"
    output += "C-term,967,1170\n"
    output += "//\n"
    output += "/path/to/Cnidaria_ginsi.fa\t# If the abs path is not defined, then the file must be in the same dir " \
              "as this .csv file.\n"
    output += "N-term,1,174\n"
    output += "TMD3,175,240\n"
    output += "ECL3,241,573\n"
    output += "TMD4,574,639\n"
    output += "ICL,640,915\n"
    output += "TMD5,916,984\n"
    output += "ECL4,985,1269\n"
    output += "TMD6,1270,1341\n"
    output += "C-term,1342,1785\n"
    output += "//\t\t# Trailing set divider is optional\n"
    print(output)
    quit()

if not in_args.blocks_file:
    sys.exit("Sorry, but you need to specify a data_blocks.csv file with -i. Use -h for instructions.")

blocks_file = os.path.abspath(in_args.blocks_file)
if not os.path.exists(blocks_file):
    sys.exit("Can't find your data_blocks.csv file at %s" % blocks_file)

with open(blocks_file, "r") as ifile:
    blocks = ifile.read()

# Clean up comments and any superfluous line breaks, then bust up the blocks into a list
blocks = re.sub("#.*", "", blocks)
blocks = re.sub("^\n", "", blocks, flags=re.MULTILINE)
blocks = re.sub("//\n", "//", blocks)
blocks = blocks.split("//")

if blocks[-1] == "":
    blocks = blocks[:-1]

csv_dir = "/".join(blocks_file.split("/")[:-1])
os.chdir(csv_dir)

# Make sure that all files specified in data_blocks actually exist, that none of those file names are duplicates, and
# that none of the sequence blocks within each set have duplicate ids.
file_names = []
for block in blocks:
    path = os.path.abspath(block.split("\n")[0])
    if not os.path.exists(path):
        sys.exit("Sorry, but you have specified a path in in your data_blocks file that does not seem to exist:\n%s  "
                 "not found." % path)

    if in_args.type == "guess":
        in_args.type = seq_tools.guess_alphabet(path)

    if in_args.type == "prot" and in_args.codon_pos:
        sys.exit("You have tried to subdevide a protein sequence up into codon positions with the -c flag.")

    file_name = path.split("/")[-1].split(".")[:-1]
    if file_name in file_names:
        sys.exit("Sorry, but you have a duplicate file name in your data_blocks file: %s" % file_name)
    else:
        file_names.append(file_name)

    seq_ranges = block.split("\n")[1:]
    if seq_ranges[-1] == "":
        seq_ranges = seq_ranges[:-1]

    range_ids = []
    for seq_range in seq_ranges:
        range_id = seq_range.split(",")[0]
        if range_id in range_ids:
            sys.exit("Sorry, but you have a duplicate data block id in the data_blocks file. The offender is under the %s set" % file_name)
        else:
            range_ids.append(range_id)

outdir = os.path.abspath(in_args.out_dir)
if not os.path.exists(outdir):
    print("Output directory not found, creating it...\n%s" % outdir)
    os.mkdir(outdir)

for block in blocks:
    path = os.path.abspath(block.split("\n")[0])
    file_name = path.split("/")[-1].split(".")[:-1]
    file_name = "_".join(file_name)
    file_name = "_".join(file_name.split(" "))

    new_dir = os.path.abspath("%s/%s" % (outdir, file_name))
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)

    infile = open(path, "r")
    phylip = open("%s/%s.phy" % (new_dir, file_name), "w")
    hashed_phylip = open("%s/%s_hashed.phy" % (new_dir, file_name), "w")
    hash_map = open("%s/%s_hash_map.csv" % (new_dir, file_name), "w")

    alignment = AlignIO.read(infile, seq_tools.guess_format(path))

    clean_alignment = seq_tools.screw_formats_align([alignment], "phylipi")
    phylip.write(clean_alignment)

    id_hashes = seq_tools.hash_seqeunce_ids(list(alignment))
    for i in id_hashes[0]:
        hash_map.write("%s,%s\n" % (i[0], i[1]))

    alignment = seq_tools.screw_formats_align([alignment], "phylipi")
    hashed_phylip.write(alignment)

    infile.close()
    phylip.close()
    hashed_phylip.close()
    hash_map.close()

    # Make cfg file
    seq_ranges = block.split("\n")[1:]
    if seq_ranges[-1] == "":
        seq_ranges = seq_ranges[:-1]

    make_cfg(new_dir, seq_ranges, file_name)
    if in_args.run_partfinder:
        os.chdir(new_dir)
        program = "partitionfinder" if in_args.type == "nucl" else "partitionfinderprotein"

        if in_args.force:
            Popen("%s --force-restart ./" % program, shell=True).wait()
        else:
            Popen("%s ./" % program, shell=True).wait()
        os.chdir(csv_dir)

print("Job complete, it ran in %s" % MyFuncs.pretty_time(round(clock()) - start_time))