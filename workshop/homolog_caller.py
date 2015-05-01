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
        self.clusters = [Cluster([y for y in x.strip().split(gene_split)]) for x in self.clusters]
        self.size = 0.
        for group in self.clusters:
            self.size += group.len

        self.printer = MyFuncs.DynamicPrint()
        self.taxa_split = taxa_split

    def score_all_clusters(self):
        score = 0
        for cluster in self.clusters:
            score += cluster.score(self.taxa_split)

        modifier = abs(self.size ** (1 / 2) - len(self.clusters))
        return round(score / len(self.clusters) - modifier ** 1.2, 2)


class Cluster():
    def __init__(self, cluster):
        self.cluster = cluster
        self.len = len(self.cluster)

    def score(self, taxa_split="-"):
        taxa = [x.split(taxa_split)[0] for x in self.cluster]
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

    def __str__(self):
        return str(self.cluster)


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

    _output = Popen("mcl %s --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" % (input_file, gq, I, tmp_dir.path),
                    shell=True, stderr=PIPE).communicate()

    _output = _output[1].decode()

    if re.search("\[mclvInflate\] warning", _output) and min_score:
        return min_score

    clusters = Clusters("%s/output.groups" % tmp_dir.path)
    score = clusters.score_all_clusters()
    return score


def homolog_caller(cluster, all_by_all, cluster_dict, rank):
    # if rank == 2:            # Delete this
    #    return cluster_dict  # Delete this

    temp_dir = MyFuncs.TempDir()

    all_by_all.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    inflation_var = mcmcmc.Variable("I", 1.4, 20)
    gq_var = mcmcmc.Variable("gq", 0.05, 0.95)

    try:
        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=100, sample_rate=1,
                                       params=["%s/input.csv" % temp_dir.path, False], outfile="%s/mcmcmc_out.csv" % temp_dir.path)
    except RuntimeError:  # Happens when mcmcmc fails to find initial chain parameters
        return cluster_dict

    # Set a 'worst score' that is reasonable for the data set
    worst_score = 10000000  # arbitrarily large number
    for chain in mcmcmc_factory.chains:
        worst_score = chain.raw_min if chain.raw_min < worst_score else worst_score

    mcmcmc_factory.reset_params(["%s/input.csv" % temp_dir.path, worst_score])

    print("\nRunning\n")
    mcmcmc_factory.run()

    mcmcmc_output = pd.read_csv("%s/mcmcmc_out.csv" % temp_dir.path, "\t")

    if rank in cluster_dict:
        cluster_dict[rank].append(cluster)
    else:
        cluster_dict[rank] = [cluster]

    best_score = max(mcmcmc_output["result"])
    if best_score < cluster.score():
        return cluster_dict

    best_df = mcmcmc_output.loc[mcmcmc_output['result'] == best_score]

    Popen("mcl %s --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" %
          ("%s/input.csv" % temp_dir.path, best_df[0:1]["gq"].values[0], best_df[0:1]["I"].values[0], temp_dir.path), shell=True, stderr=PIPE).communicate()

    mcl_clusters = Clusters("%s/output.groups" % temp_dir.path)
    for clust in mcl_clusters.clusters:
        if len(clust.cluster) in [1, 2]:
            if rank + 1 in cluster_dict:
                cluster_dict[rank + 1].append(clust)
            else:
                cluster_dict[rank + 1] = [clust]
            continue

        group_all_by_all = split_all_by_all(all_by_all, clust.cluster)["removed"]

        _min = group_all_by_all[:][2].min()
        data_range = group_all_by_all[:][2].max() - _min

        re_norm = (group_all_by_all[:][2] - _min) / data_range

        group_all_by_all[:][2] = re_norm

        cluster_dict = homolog_caller(clust, group_all_by_all, cluster_dict, rank + 1)

    return cluster_dict

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="homolog_caller", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("all_by_all_file", help="Location of all-by-all scores file", action="store")
    parser.add_argument("-sep", "--separator", action="store", default="\t",
                        help="If the all-by-all file is not tab delimited, specify the character")

    in_args = parser.parse_args()

    best = None
    best_clusters = None
    lock = Lock()

    counter = 1

    scores_data = pd.read_csv(os.path.abspath(in_args.all_by_all_file), in_args.separator, header=None)

    master_cluster = pd.concat([scores_data[0], scores_data[1]])
    master_cluster = master_cluster.value_counts()
    master_cluster = Cluster([i for i in master_cluster.index])

    final_cluster_dict = {}
    final_cluster_dict = homolog_caller(master_cluster, scores_data, final_cluster_dict, 0)

    output = ""
    for i in final_cluster_dict:
        for j in final_cluster_dict[i]:
            for k in j.cluster:
                output += "%s\t" % k
            output = "%s\n" % output.strip()
        output += "\n"

    with open("testoutput.txt", "w") as ofile:
        ofile.write(output)