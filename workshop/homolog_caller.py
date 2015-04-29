#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 27 2015 

"""
DESCRIPTION OF PROGRAM
"""
# Std library
import os
import re
from copy import copy
from subprocess import Popen, PIPE
from multiprocessing import Lock

# 3rd party
import pandas as pd

# My packages
import mcmcmc
import MyFuncs


class Clusters():
    def __init__(self, path, group_split="\n", gene_split="\t", taxa_split="-"):
        with open(path, "r") as ifile:
            self.input = ifile.read()

        self.clusters = self.input.strip().split(group_split)
        self.clusters = [[y for y in x.strip().split(gene_split)] for x in self.clusters]
        self.clusters.reverse()
        self.size = 0.
        for group in self.clusters:
            self.size += len(group)
        self.printer = MyFuncs.DynamicPrint()
        self.taxa_split = taxa_split

    def score_all_clusters(self):
        score = 0
        for cluster in self.clusters:
            score += self.score_cluster(cluster, self.taxa_split)

        modifier = abs(self.size ** (1 / 2) - len(self.clusters))
        return round(score / len(self.clusters) - modifier ** 1.2, 2)

    @staticmethod
    def score_cluster(cluster, taxa_split="-"):
        taxa = [x.split(taxa_split)[0] for x in cluster]
        taxa = pd.Series(taxa)

        if len(taxa) == 1:
            return 0

        score = 0
        singles = 0
        for taxon in taxa.value_counts():
            if taxon == 1:
                singles += 1

            elif taxon > 1:
                score -= taxon

            else:
                raise ValueError("A taxon was found with %s records" % taxon)

        score += singles ** 2
        return score


def split_all_by_all(data_frame, remove_list):
    if len(data_frame.columns) != 3:
        raise AttributeError("dataframe should be 3 columns")

    data_frame = copy(data_frame)
    columns = data_frame.columns
    data_frame.columns = [0, 1, 2]

    removed = data_frame.loc[data_frame[0].isin(remove_list)]
    removed = removed.loc[removed[1].isin(remove_list)]

    remaining = data_frame.loc[~data_frame[0].isin(remove_list)]
    remaining = remaining.loc[~remaining[1].isin(remove_list)]

    remaining.columns = columns
    removed.columns = columns
    return {"remaining": remaining, "removed": removed}


def mcmcmc_mcl(args, params):
    I, gq = args
    input_file, min_score = params
    tmp_dir = MyFuncs.TempDir()

    output = Popen("mcl %s --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" % (input_file, gq, I, tmp_dir.path),
                   shell=True, stderr=PIPE).communicate()

    output = output[1].decode()

    if re.search("\[mclvInflate\] warning", output) and min_score:
        return min_score

    clusters = Clusters("%s/output.groups" % tmp_dir.path)
    score = clusters.score_all_clusters()
    return score


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="homolog_caller", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("all_by_all_file", help="Location of all-by-all scores file", action="store")
    parser.add_argument("-sep", "--separator", action="store", default="\t",
                        help="If the all-by-all file is not tab delimited, specify the character")

    in_args = parser.parse_args()

    scores_data = pd.read_csv(os.path.abspath(in_args.all_by_all_file), in_args.separator, header=None)

    best = None
    best_clusters = None
    lock = Lock()
    inflation_steps = [20, 10, 6, 5, 4, 3, 2, 1.4]
    gq_vals = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9]
    mcl_params = []
    temp_dir = MyFuncs.TempDir()
    counter = 1

    for inflation in inflation_steps:
        for gq in gq_vals:
            mcl_params.append({"work_dir": temp_dir.path, "I": inflation,
                               "gq": "'gq(%s)'" % gq, "out_name": "mcl%s" % counter})
            counter += 1

    scores_data.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    inflation_var = mcmcmc.Variable("I", 1.4, 20)
    gq_var = mcmcmc.Variable("gq", 0.05, 0.95)

    mcmcmc = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=10, sample_rate=1,
                           params=["%s/input.csv" % temp_dir.path, False])

    # Set a 'worst score' that is reasonable for the data set
    worst_score = 10000000  # arbitrarily large number
    for chain in mcmcmc.chains:
        worst_score = chain.raw_min if chain.raw_min < worst_score else worst_score

    mcmcmc.reset_params(["%s/input.csv" % temp_dir.path, worst_score])
    print(worst_score)
    print("\nRunning\n")
    mcmcmc.run()
