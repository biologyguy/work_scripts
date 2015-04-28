#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 27 2015 

"""
DESCRIPTION OF PROGRAM
"""

import sys
import os
import re
#import shutil
import MyFuncs
#import SeqBuddy
import pandas as pd
from copy import copy
from subprocess import Popen, PIPE
from multiprocessing import Lock


def fibonacci(val):
    lookup = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711,
              28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887]
    try:
        return lookup[val]
    except IndexError:
        return fibonacci(val - 1) + fibonacci(val - 2)


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
        #print(taxa.value_counts())
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


def run_mcl(scores_df):
    temp_dir = MyFuncs.TempDir()
    scores_df.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    clusters = Clusters("%s/mclOutput" % temp_dir.path)

    return


def mc_mcl(args):
    if not args["gq"]:
        mcl_output = Popen("mcl %s/input.csv --abc -te 2 -I %s -o %s/%s.groups" %
                           (args["work_dir"], args["I"], args["work_dir"], args["out_name"]),
                           shell=True, stderr=PIPE).communicate()

    else:
        mcl_output = Popen("mcl %s/input.csv --abc -te 2 -tf %s -I %s -o %s/%s.groups" %
                           (args["work_dir"], args["gq"], args["I"], args["work_dir"], args["out_name"]),
                           shell=True, stderr=PIPE).communicate()

    with open("%s/%s.stdout" % (args["work_dir"], args["out_name"]), "w") as _ofile:
        _ofile.write(mcl_output[1].decode())


    output_txt = mcl_output[1].decode()

    if re.search("\[mclvInflate\] warning", output_txt):
        return

    num_clusters = int(re.search("\[mcl\] ([0-9]+) clusters found", output_txt).group(1))

    if num_clusters == 1:
        return

    clusters = Clusters("%s/%s.groups" % (args["work_dir"], args["out_name"]))
    score = clusters.score_all_clusters()
    args['gq'] = re.search("gq\((0\.[0-9]+)\)", args['gq']).group(1)
    with lock:
        print("%s\t%s\t%s\t%s" % (args["I"], round(float(args["gq"]), 2), num_clusters, "{:,}".format(score)))
    return



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="homolog_caller", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("all_by_all_file", help="Location of all-by-all scores file", action="store")
    parser.add_argument("-sep", "--separator", action="store", default="\t",
                        help="If the all-by-all file is not tab delimited, specify the character")
    #parser.add_argument("-c", "--choice", help="", type=str, choices=["", ""], default=False)
    #parser.add_argument("-m", "--multi_arg", nargs="+", help="", default=[])

    in_args = parser.parse_args()

    scores_data = pd.read_csv(os.path.abspath(in_args.all_by_all_file), in_args.separator, header=None)

    """
    Run through 4 mc_mcl() sets:
        1) Check -I between 10 and 100
        2) Check -I between the best set from previous step
        3) Check -gp between 0.05 and 0.95 with best -I
        4) Check -gp between the best set form previous step
    """

    best = None
    best_clusters = None
    print("Inf\tgq\t#groups\tscore")
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
    MyFuncs.run_multicore_function(mcl_params, mc_mcl, quiet=True)

    """
    for i in range(0, 10):
        with open("%s/mcl%s.stdout" % (temp_dir.path, i), "r") as ifile:
            output_txt = ifile.read()

        if re.search("\[mclvInflate\] warning", output_txt):
            continue

        num_clusters = int(re.search("\[mcl\] ([0-9]+) clusters found", output_txt).group(1))

        if num_clusters == 1:
            continue

        clusters = Clusters("%s/mcl%s.groups" % (temp_dir.path, i))
        score = clusters.score_all_clusters()
        print("%s\t%s\t%s" % (inflation, gq_vals[i], "{:,}".format(score)))


        score = 1
        for clust in clusters.clusters:
            score *= len(clust)

        max_score = (clusters.size ** (1 / 2)) ** (clusters.size ** (1 / 2))
        score = int(abs(max_score - score))
        print("%s\t%s\t%s" % (inflation, gq_vals[i], score))

        if not best or max_score - score < best:
            best = num_clusters
            best_clusters = Clusters("%s/mcl%s.groups" % (temp_dir.path, i))
        """
    #clusters = Clusters("%s/mclOutput" % temp_dir.path)
