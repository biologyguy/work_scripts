#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Feb 23 2015 

"""
DESCRIPTION OF PROGRAM
"""

#import sys
import re
import argparse
import SeqBuddy

parser = argparse.ArgumentParser(prog="test_panx_file", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("sequences", help="Sequence file to compare", action="store")
parser.add_argument("id_file", help="File containing all sequence ids", action="store")

in_args = parser.parse_args()

if __name__ == '__main__':
    sequences = SeqBuddy.SeqBuddy(in_args.sequences)
    ids_in_seq_file = [rec.id for rec in sequences.records]
    panx_list = []

    with open(in_args.id_file, "r") as ifile:
        ids_file = ifile.read().split("\n\n")

    print("IDs in name map, not in sequence file.")
    for taxa in ids_file:
        lines = taxa.split("\n")
        species = re.search('- (.*)', lines[0]).group(1)
        printed = False
        for panx in lines[1:]:
            panx = panx.split("\t")[-1]
            panx_list.append(panx)
            if panx not in ids_in_seq_file:
                if not printed:
                    printed = True
                    print(species)
                print(panx)

        if printed:
            print("")

    print("\nIDs in sequence file, not in name map.")
    for panx in ids_in_seq_file:
        if panx not in panx_list:
            print(panx)
