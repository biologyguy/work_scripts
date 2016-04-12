#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Mar 21 16

import SeqBuddy as Sb
import AlignBuddy as Alb
from Bio.SubsMat import SeqMat, MatrixInfo
from Bio.Seq import Seq
from copy import copy
from math import log
from random import sample, randint
import re
import MyFuncs
from multiprocessing import Lock
import sys
import os
import numpy as np


def make_full_mat(subsmat):
    for key in copy(subsmat):
        try:
            # don't over-write the reverse keys if they are already initialized
            subsmat[(key[1], key[0])]
        except KeyError:
            subsmat[(key[1], key[0])] = subsmat[key]
    return subsmat


def bit_score(raw_score):
    # These values were empirically determined for BLOSUM62 by Altschul
    bit_k_value = 0.035
    bit_lambda = 0.252

    bits = ((bit_lambda * raw_score) - (log(bit_k_value))) / log(2)
    return bits


def shuffle(_seq):
    _tmds = ""
    _other = ""
    for _feat in _seq.features:
        _extract = _feat.extract(_seq)
        if "TMD" in _feat.type:
            _tmds += str(_extract.seq)
        else:
            _other += str(_extract.seq)

    _tmds = ''.join(sample(_tmds, len(_tmds)))
    _other = ''.join(sample(_other, len(_other)))

    new_seq = ""
    for _feat in _seq.features:
        _feat_len = len(_feat.extract(_seq))
        if "TMD" in _feat.type:
            new_seq += _tmds[:_feat_len]
            _tmds = _tmds[_feat_len:]
        else:
            new_seq += _other[:_feat_len]
            _other = _other[_feat_len:]

    _seq.seq = Seq(new_seq, _seq.seq.alphabet)
    return _seq


def score_sequences(seq_pair):
    seq1, seq2 = seq_pair.records
    id_regex = "^%s$|^%s$" % (seq1.id, seq2.id)
    sb_copy = Sb.make_copy(seqbuddy)
    Sb.delete_records(sb_copy, id_regex)
    sb_copy = Sb.SeqBuddy(sb_copy.records + [seq1, seq2], out_format="gb", alpha=sb_copy.alpha)
    alignbuddy = Alb.generate_msa(sb_copy, tool="mafft", params=" --globalpair", quiet=True)
    if not in_args.no_msa_trim:
        alignbuddy = Alb.trimal(alignbuddy, threshold="gappyout")
    alignbuddy = Alb.pull_records(alignbuddy, id_regex)
    _score = 0
    seq1, seq2 = alignbuddy.records()
    prev_aa1 = "-"
    prev_aa2 = "-"

    for aa_pos in range(alignbuddy.lengths()[0]):
        aa1 = seq1.seq[aa_pos]
        aa2 = seq2.seq[aa_pos]

        if aa1 == "-" or aa2 == "-":
            if prev_aa1 == "-" or prev_aa2 == "-":
                _score += gap_extend
            else:
                _score += gap_open
        else:
            _score += BLOSUM45[aa1, aa2]
        prev_aa1 = str(aa1)
        prev_aa2 = str(aa2)
    return _score


def mc_run_comparison(seq_pair, args):
    # pairwise alignments

    temp_file_path = args[0]
    seq_pair_sb = Sb.pull_recs(Sb.make_copy(seqbuddy), "^%s$|^%s$" % (seq_pair[0], seq_pair[1]))

    score = score_sequences(seq_pair_sb)
    shuffle_scores = []
    for _ in range(in_args.num_shuffles):
        seq_pair_copy = Sb.make_copy(seq_pair_sb)
        seq1 = shuffle(seq_pair_copy.records[0])
        seq2 = shuffle(seq_pair_copy.records[1])
        shuffle_scores.append(score_sequences(Sb.SeqBuddy([seq1, seq2])))

    ave_shuffles = np.average(shuffle_scores)
    stdev_shuffles = np.std(shuffle_scores)
    with lock:
        with open(temp_file_path, "a") as ofile:
            ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seq_pair[0], seq_pair[1], score, ave_shuffles,
                                                              stdev_shuffles, score - ave_shuffles, bit_score(score),
                                                              bit_score(ave_shuffles)))
    return

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="ancestral_recon_compare.py", description="")

    parser.add_argument("in_file", help="Location of genbank input file", action="store")
    parser.add_argument("-n", "--num_shuffles", help="How many shuffles to do", type=int, default=10)
    parser.add_argument("-op", "--open_penalty", help="Penalty for opening a gap in pairwise alignment scoring",
                        type=float, default=-10)
    parser.add_argument("-ep", "--extend_penalty", help="Penalty for extending a gap in pairwise alignment scoring",
                        type=float, default=-1)
    parser.add_argument("-nt", "--no_msa_trim", action="store_true",
                        help="Don't apply the gappyout algorithm to MSAs before scoring")

    in_args = parser.parse_args()

    PHAT = make_full_mat(SeqMat(MatrixInfo.phat75_73))
    BLOSUM62 = make_full_mat(SeqMat(MatrixInfo.blosum62))
    BLOSUM45 = make_full_mat(SeqMat(MatrixInfo.blosum45))

    gap_open = in_args.open_penalty if in_args.open_penalty <= 0 else in_args.open_penalty * -1
    gap_extend = in_args.extend_penalty if in_args.extend_penalty <= 0 else in_args.extend_penalty * -1
    shuffle_iters = in_args.num_shuffles
    input_file = os.path.abspath(in_args.in_file)

    temp_file = MyFuncs.TempFile()
    lock = Lock()
    # alignbuddy = Alb.AlignBuddy(input_file)
    seqbuddy = Sb.SeqBuddy(input_file)

    # Need to add annotations for non-TMD regions (used later)
    for rec in seqbuddy.records:
        new_feats = []
        previous_start = 1
        for feat in rec.features:
            new_feats.append("%s-%s" % (previous_start, feat.location.start))
            previous_start = feat.location.end + 1
        new_feats.append("%s-%s" % (previous_start, len(rec.seq)))
        temp_seqbuddy = Sb.SeqBuddy([rec])
        for feat in new_feats:
            Sb.annotate(temp_seqbuddy, _type="other", location=feat)

    anc_panx_ids = []
    panx_ids = []
    anc_cx_ids = []
    cx_ids = []

    for rec in seqbuddy.records:
        if re.match("Anc_Panx", rec.id):
            anc_panx_ids.append(rec.id)
        elif re.match("Anc_Cx", rec.id):
            anc_cx_ids.append(rec.id)
        elif re.search("\-Panx", rec.id):
            panx_ids.append(rec.id)
        elif re.search("\-Cx", rec.id):
            cx_ids.append(rec.id)
        else:
            raise IndexError("Not sure what to do with record ID '%s'" % rec.id)

    sequence_pairs = []
    panx_copy = list(anc_panx_ids + panx_ids)
    for panx1 in anc_panx_ids + panx_ids:
        panx_copy.remove(panx1)
        for panx2 in panx_copy:
            sequence_pairs.append((panx1, panx2))
        for cx in anc_cx_ids + cx_ids:
            sequence_pairs.append((panx1, cx))

    MyFuncs.run_multicore_function(sequence_pairs, mc_run_comparison, [temp_file.path], out_type=sys.stderr)

    with open(temp_file.path, "r") as ifile:
        panxs_to_cxs = ifile.read().split("\n")

    print("Panx\tCx\tScore\tShuffle_ave\tShuffle_stdev\tDelta\tBit_score\tShuff_bit_score")

    for line in panxs_to_cxs:
        print(line)
