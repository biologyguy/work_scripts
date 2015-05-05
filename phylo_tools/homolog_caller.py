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
        self.cluster = cluster.sort()
        self.len = len(self.cluster)
        self.name = ""

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


def homolog_caller(cluster, all_by_all, cluster_list, rank, save=False, steps=1000):
    temp_dir = MyFuncs.TempDir()

    all_by_all.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    inflation_var = mcmcmc.Variable("I", 1.4, 20)
    gq_var = mcmcmc.Variable("gq", 0.05, 0.95)

    try:
        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1,
                                       params=["%s/input.csv" % temp_dir.path, False], quiet=True,
                                       outfile="%s/mcmcmc_out.csv" % temp_dir.path)
    except RuntimeError:  # Happens when mcmcmc fails to find different initial chain parameters
        cluster_list.append(cluster)
        if save:
            temp_dir.save("%s/group_%s" % (save, rank))
        return cluster_list

    # Set a 'worst score' that is reasonable for the data set
    worst_score = 10000000  # arbitrarily large number to start
    for chain in mcmcmc_factory.chains:
        worst_score = chain.raw_min if chain.raw_min < worst_score else worst_score

    mcmcmc_factory.reset_params(["%s/input.csv" % temp_dir.path, worst_score])

    mcmcmc_factory.run()

    mcmcmc_output = pd.read_csv("%s/mcmcmc_out.csv" % temp_dir.path, "\t")

    best_score = max(mcmcmc_output["result"])
    if best_score < cluster.score():
        cluster_list.append(cluster)
        if save:
            temp_dir.save("%s/group_%s" % (save, rank))
        return cluster_list

    best_df = mcmcmc_output.loc[mcmcmc_output['result'] == best_score]

    Popen("mcl %s --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" %
          ("%s/input.csv" % temp_dir.path, best_df[0:1]["gq"].values[0], best_df[0:1]["I"].values[0], temp_dir.path), shell=True, stderr=PIPE).communicate()

    mcl_clusters = Clusters("%s/output.groups" % temp_dir.path)
    _counter = 1
    for clust in mcl_clusters.clusters:
        next_rank = "%s_%s" % (rank, _counter)
        clust.name = next_rank
        _counter += 1
        if len(clust.cluster) in [1, 2]:
            cluster_list.append(clust)
            continue

        group_all_by_all = split_all_by_all(all_by_all, clust.cluster)["removed"]

        _min = group_all_by_all[:][2].min()
        data_range = group_all_by_all[:][2].max() - _min

        re_norm = (group_all_by_all[:][2] - _min) / data_range

        group_all_by_all[:][2] = re_norm

        # Recursion...
        cluster_list = homolog_caller(clust, group_all_by_all, cluster_list, next_rank, save, steps=steps)

    if save:
        temp_dir.save("%s/group_%s" % (save, rank))

    return cluster_list

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="homolog_caller", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("all_by_all_file", help="Location of all-by-all scores file", action="store")
    parser.add_argument("output_file", help="Where should groups be written to? Default to stdout.", action="store", nargs="?")
    parser.add_argument("-mcs", "--mcmcmc_steps", default=1000, type=int,
                        help="Specify how deeply to sample MCL parameters")
    parser.add_argument("-sep", "--separator", action="store", default="\t",
                        help="If the all-by-all file is not tab delimited, specify the character")
    parser.add_argument("-stf", "--save_temp_files", help="Keep all mcmcmc and MCL files", default=False)

    in_args = parser.parse_args()

    best = None
    best_clusters = None
    lock = Lock()

    counter = 1

    scores_data = pd.read_csv(os.path.abspath(in_args.all_by_all_file), in_args.separator, header=None)

    master_cluster = pd.concat([scores_data[0], scores_data[1]])
    master_cluster = master_cluster.value_counts()
    master_cluster = Cluster([i for i in master_cluster.index])

    if in_args.save_temp_files:
        if not os.path.isdir(in_args.save_temp_files):
            os.makedirs(in_args.save_temp_files)

    print("Executing Homolog Caller...")
    final_clusters = []
    final_clusters = homolog_caller(master_cluster, scores_data, final_clusters, 0,
                                    in_args.save_temp_files, steps=in_args.mcmcmc_steps)

    output = ""
    for clust in final_clusters:
        output += "group_%s\t%s\t" % (clust.name, clust.score())
        for seq_id in clust.cluster:
            output += "%s\t" % seq_id
        output = "%s\n" % output.strip()

    if not in_args.output_file:
        print(output)

    else:
        with open(in_args.output_file, "w") as ofile:
            ofile.write(output)