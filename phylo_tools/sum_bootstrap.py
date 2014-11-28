#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 28 2014 

"""
Silly little script to get some data on the bootstrap score of a tree
"""

import re
import sys
import argparse

parser = argparse.ArgumentParser(prog="sum_bootstrap", description="Add up all the bootstrap scores in a tree file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("newick_file", help="What tree do you want analysed?", action="store")

in_args = parser.parse_args()


with open(in_args.newick_file, "r") as ifile:
    content = ifile.read()

# this is for 'raw' newick, as output by RAxML
bootstraps = re.finditer("\)([0-9]{1,3}):", content)
boot_sum = 0
counter = 0
for boot in bootstraps:
    boot_sum += int(boot.group(1))
    counter += 1

# If the newick has been modified by FigTree...
if counter == 0:
    bootstraps = re.finditer("bootstrap[s]*=([0-9]{1,3})", content)
    for boot in bootstraps:
        boot_sum += int(boot.group(1))
        counter += 1

if counter == 0:
    print("Unable to find the bootstrap values in your file.")
else:
    print("Sum:\t%s\nAve:\t%s" % (boot_sum, round(boot_sum / float(counter), 2)))