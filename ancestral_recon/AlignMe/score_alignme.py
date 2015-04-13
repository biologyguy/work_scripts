#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re


def score_alignme(alignme_file):
    with open(alignme_file, "r") as ifile:
        file_lines = ifile.readlines()

    # clear out header rows
    while True:
        if not re.match("#", file_lines[0]):
            del file_lines[0]
        else:
            break
    del file_lines[-1]

    # determine normalizing maximum. I have semi-arbitrarily set this as one-half the sum of residues in both seqs
    seq1_len = 0.0
    seq2_len = 0.0
    for data in file_lines:
        data = re.split("\s+", data)
        if data[1] != "?0":
            seq1_len += 1
        if data[5] != "?0":
            seq2_len += 1

    normalizing_len = (seq1_len + seq2_len) / 2

    tally = 0.0
    for data in file_lines:
        if data[0] == "#":
            continue
        data = re.sub("\s+", ",", data)
        data = data.split(",")
        if "?0" in [data[1], data[5]]:  # Gaps get a score of 0
            continue

        # At the moment, all four parameters are given equal weight
        membrane_score = 1 - abs(float(data[1]) - float(data[5]))
        coil_score = 1 - abs(float(data[2]) - float(data[6]))
        helix_score = 1 - abs(float(data[3]) - float(data[7]))
        sheet_score = 1 - abs(float(data[4]) - float(data[8]))
        tally += membrane_score + coil_score + helix_score + sheet_score

    return round(tally / normalizing_len, 5) / 4  # The '4' is for the number of columns being compared

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="", description="")
    parser.add_argument('in_file', help='', action='store')
    in_args = parser.parse_args()

    print(score_alignme(in_args.in_file))