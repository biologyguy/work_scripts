#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 22 2015 

"""
DESCRIPTION OF PROGRAM
"""

#import sys
#import os
#import re
#import shutil
import MyFuncs
#import SeqBuddy
#import argparse


class Clusters():
    def __init__(self, path):
        with open(path, "r") as ifile:
            self.input = ifile.read()

        self.clusters = self.input.strip().split("group")[1:]
        self.clusters = [[y for y in x.strip().split(" ")[1:]] for x in self.clusters]
        self.size = 0.
        for group in self.clusters:
            self.size += len(group)
        self.printer = MyFuncs.DynamicPrint()

    def compare(self, query_clusters):
        score = 0.
        counter = 1
        for subj in self.clusters:
            printer.write("Cluster %s of %s" % (counter, len(self.clusters)))
            counter += 1
            tally = 0.
            len_subj = len(subj)

            for query in query_clusters.clusters:
                matches = self.num_matches(subj, query)
                if not matches:
                    continue
                else:
                    tally += (matches * 2.) / (len(subj) + len(query))
                    len_subj -= matches
                    if len_subj == 0:
                        score += tally * (len(subj) / self.size)
                        break
        print("")
        return score

    @staticmethod
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

    timer = MyFuncs.Timer()
    printer = MyFuncs.DynamicPrint()


    groups1 = Clusters("test_files/e_10_groups.txt.bak")
    groups2 = Clusters("test_files/e_10-5_groups.txt.bak")
    groups3 = Clusters("test_files/test_groups3.txt")

    print("Score: %s\n%s" % (groups1.compare(groups2), timer.end()))

    """
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
    """
