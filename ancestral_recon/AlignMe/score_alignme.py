#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys


def score_alignme(alignme_file):
    file_list = alignme_file.readlines()

    # clear out header rows
    while True:
        if not re.match("#", file_list[0]):
            del file_list[0]
        else:
            break
    del file_list[-1]

    # determine normalizing maximum
    seq1_len = 0.0
    seq2_len = 0.0
    for _next in file_list:
        data = re.split("\s+", _next)
        if data[1] != "?0":
            seq1_len += 1
        if data[5] != "?0":
            seq2_len += 1

    norm_max = (seq1_len + seq2_len) * 4
    print(norm_max)

    tally = 0.0
    count = 0
    for _next in file_list:
        regular = re.sub("\s+", ",", _next)
        regular = re.sub("\?0", "0", regular)
        data = regular.split(",")
        tally += abs(float(data[1]) - float(data[5])) + abs(float(data[2]) - float(data[6])) + \
                 abs(float(data[3]) - float(data[7])) + abs(float(data[4]) - float(data[8]))

    return round(tally / norm_max, 3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="", description="")
    parser.add_argument('in_file', help='', action='store')
    in_args = parser.parse_args()

    with open(in_args.in_file, "r") as ifile:
        print(score_alignme(ifile))