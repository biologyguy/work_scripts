#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 22 2015 

"""
DESCRIPTION OF PROGRAM
"""

import sys
import os
import re
import shutil
import MyFuncs
import SeqBuddy
import argparse


class NewClass():
    """DESCRIPTION OF CLASS"""
    def __init__(self):
        self.x = 1

    def class_def(self):
        self.x = 1
        return self.x


def parse_clusters(file_handle):
    """
    Create a dictionary with indices for every sequence, containing all other sequences in their group
    """
    clusters = {}
    for group in file_handle:
        genes = group.split("\t")
        for gene in genes:
            clusters[gene.strip()] = [x.strip() for x in genes]
    return clusters


def num_matches(_subj, _query):
    count = 0.
    for next_item in _subj:
        if next_item in _query:
            count += 1
    return count if count > 0 else None


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(prog="compare_homolog_groups", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument("positional_arg1", help="", action="store")
    # parser.add_argument("-t", "--true", help="", action="store_true", default=False)
    # parser.add_argument("-c", "--choice", help="", type=str, choices=["", ""], default=False)
    # parser.add_argument("-m", "--multi_arg", nargs="+", help="", default=[])
    # in_args = parser.parse_args()

    printer = MyFuncs.DynamicPrint()

    with open("test_files/e_10_groups.txt", "r") as ifile:
        groups1 = ifile.read()
        groups1 = groups1.split("group")
        groups1 = [[y for y in x.strip().split(" ")[1:]] for x in groups1]

    with open("test_files/e_10-5_groups.txt", "r") as ifile:
        groups2 = ifile.read()
        groups2 = groups2.split("group")
        groups2 = [[y for y in x.strip().split(" ")[1:]] for x in groups2]

    g1_size = 0.
    for group in groups1:
        g1_size += len(group)

    g2_size = 0.
    for group in groups2:
        g2_size += len(group)

    g1_score = 0.
    for subj in groups1:
        best = 0.
        len_subj = len(subj)

        for query in groups2:
            matches = num_matches(subj, query)
            if not matches:
                continue
            else:
                score = (matches * 2.) / (len(subj) + len(query))
                best = score if score > best else best
                len_subj -= matches

                if len_subj == 0:
                    g1_score += best / g1_size
                    break

    g2_score = 0.
    for subj in groups2:
        best = 0.
        len_subj = len(subj)

        for query in groups1:
            matches = num_matches(subj, query)
            if not matches:
                continue
            else:
                score = (matches * 2.) / (len(subj) + len(query))
                best = score if score > best else best
                len_subj -= matches

                if len_subj == 0:
                    g2_score += best / g1_size
                    break

    print("G1: %s, G2: %s" % (g1_score, g2_score))