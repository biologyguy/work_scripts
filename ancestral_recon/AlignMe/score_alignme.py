#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from math import log
from SeqBuddy import SeqBuddy
import sys
from Bio.SubsMat import MatrixInfo, SeqMat


PHAT = SeqMat(MatrixInfo.phat75_73)
BLOSUM62 = SeqMat(MatrixInfo.blosum62)

def bit_score(raw_score):
    # These values were empirically determined for BLOSUM62 by Altschul
    bit_k_value = 0.035
    bit_lambda = 0.252

    bits = ((bit_lambda * raw_score) - (log(bit_k_value))) / log(2)
    return bits


def scale_range(data, percentile=1.0):
    if 0. > percentile > 1.0:
        raise ValueError("scale_range() percentile parameter should be between 0.5 and 1.0")

    if percentile < 0.5:
        percentile = 1 - percentile

    max_limit = round(len(data) * percentile) - 1
    min_limit = -1 * (max_limit + 1)

    if type(data) == dict:
        sorted_data = sorted([data[key] for key in data])
        _max = sorted_data[max_limit]
        _min = sorted_data[min_limit]
        data_range = _max - _min
        for key in data:
            data[key] = (data[key] - _min) / data_range
            data[key] = 1. if data[key] > 1. else data[key]
            data[key] = 0. if data[key] < 0. else data[key]

    else:
        sorted_data = sorted(data)
        _max = sorted_data[max_limit]
        _min = sorted_data[min_limit]
        data_range = _max - _min
        for i in range(len(data)):
            data[i] = (data[i] - _min) / data_range
            data[i] = 1. if data[i] > 1. else data[i]
            data[i] = 0. if data[i] < 0. else data[i]

    return data


def alignment_sub_mat_score(subj_top, query_top, subj_align, query_align):
    normalizing_len = (len(subj_top) + len(query_top)) / 2
    score = 0
    gaps = 0
    gaps = 0
    for i in range(len(subj_align)):
        if subj_align[i] == "-":
            gaps += 1
            query_top = query_top[1:]

        elif query_align[i] == "-":
            gaps += 1
            subj_top = subj_top[1:]

        else:
            _pair = sorted((subj_align[i], query_align[i]))
            _pair = tuple(_pair)
            if subj_top[0] == "M" and query_top[0] == "M":
                score += PHAT[_pair]
            else:
                score += BLOSUM62[_pair]

            subj_top = subj_top[1:]
            query_top = query_top[1:]

    score = bit_score(score)
    score /= normalizing_len
    return {"score": score, "gaps": gaps}

if __name__ == '__main__':

    with open("/Users/bondsr/Documents/work_scripts/ancestral_recon/AlignMe/Testing/ALIGNME_FILES/Cel-Panxε6-Mle-Panxα9.aln", 'r') as ifile:
        data = ifile.read()

    regex = re.sub("^[^ ].* +$", '', data, flags=re.MULTILINE)
    print(regex)
    #print(data)
    sys.exit()
    tops = SeqBuddy("/Users/bondsr/Documents/work_scripts/ancestral_recon/AlignMe/Alignme_output/AligneMe_input_top.fasta")

    pairwise_array = []
    for x in tops.records:
        for y in tops.records:
            if (x.id, y.id) in pairwise_array or (y.id, x.id) in pairwise_array or x.id == y.id:
                continue
            else:
                pairwise_array.append((x.id, y.id))

    tops = tops.to_dict()
    scores = {}
    for pair in pairwise_array:
        align = SeqBuddy("/Users/bondsr/Documents/work_scripts/ancestral_recon/AlignMe/Alignme_output/ALIGNME_FILES/%s-%s.aln" % (pair[0], pair[1]), "clustal")
        align = align.records
        subj_align = str(align[0].seq)
        query_align = str(align[1].seq)

        subj_top = str(tops[pair[0]].seq)
        query_top = str(tops[pair[1]].seq)

        output = alignment_sub_mat_score(subj_top, query_top, subj_align, query_align)
        scores["%s-%s" % (pair[0], pair[1])] = output["score"]
        # print("%s-%s\t%s\t%s" % (pair[0], pair[1], round(output["score"], 4), output["gaps"]))

    scores = scale_range(scores, 0.95)
    for key in scores:
        print("%s:\t%s" % (key, round(scores[key], 3)))
    import numpy as np


    # for pair in scores:
    #     print("%s\t%s" % (pair, round(scores[pair], 4)))


"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="", description="")
    parser.add_argument('in_file', help='', action='store')
    in_args = parser.parse_args()

    print(score_alignme(in_args.in_file))
"""