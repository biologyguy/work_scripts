#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Oct 24 2014 

"""
DESCRIPTION OF PROGRAM
"""

import argparse
import re, os, sys

parser = argparse.ArgumentParser(prog="best_ML_tree", description="Glean the best RAxML tree from the RAxML_info file after a run")

parser.add_argument("info_file", help="Where is the info file?", action="store")
parser.add_argument("-f", "--force", help="Overwrite RAxML_best. file if exists. Careful...", default=False)

in_args = parser.parse_args()
info_file = os.path.abspath(in_args.info_file)

with open(info_file, "r") as ifile:
    content = ifile.read()

inferences = re.findall("Inference\[[0-9]*\] final.*", content)

best_score = -9999999999999
best_index = 0
tree_file = ""
for i in range(len(inferences)):
    line = inferences[i]
    if line == "":
        break
    chunks = line.split(" ")
    if float(chunks[4]) > best_score:
        best_score = float(chunks[4])
        best_index = i
        tree_file = chunks[9]

if not os.path.exists(tree_file):
    sys.exit("ERROR: Could not find tree file %s" % tree_file)

info_file_path = info_file.split("/")
outpath = "/".join(info_file_path[:-1])
extension = info_file_path[-1].split(".")[-1]
outfile = "%s/RAxML_best.%s" % (outpath, extension)

if os.path.exists(outfile) and not in_args.force:
    sys.exit("WARNING: %s already exists, use -f flag to overwrite" % outfile)

with open(tree_file, "r") as ifile:
    content = ifile.read()

with open(outfile, "w") as ofile:
    ofile.write(content)