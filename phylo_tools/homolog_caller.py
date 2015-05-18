#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 27 2015 

"""
DESCRIPTION OF PROGRAM
"""
# Std library
import sys
import os
import re
from copy import copy
from subprocess import Popen, PIPE
from multiprocessing import Lock
from random import sample
from math import floor

# 3rd party
import pandas as pd
import numpy as np
import statsmodels.api as sm

# My packages
import mcmcmc
import MyFuncs


class Clusters():
    def __init__(self, path, group_split="\n", gene_split="\t", taxa_split="-"):
        with open(path, "r") as _ifile:
            self.input = _ifile.read()

        self.clusters = self.input.strip().split(group_split)
        self.clusters = [Cluster([y for y in x.strip().split(gene_split)]) for x in self.clusters]
        self.size = 0.
        for group in self.clusters:
            self.size += group.len

        self.printer = MyFuncs.DynamicPrint()
        self.taxa_split = taxa_split

    def score_all_clusters(self):
        score = 0
        taxa = pd.Series()
        for cluster in self.clusters:
            # Weight each score by cluster size. Final score max = 1.0, min = 0.0
            score += cluster.score(self.taxa_split) * (cluster.len / self.size)
            temp_taxa = [x.split(self.taxa_split)[0] for x in cluster.cluster]
            temp_taxa = pd.Series(temp_taxa)
            taxa = pd.concat([taxa, temp_taxa])

        # print(modifier)
        # print(taxa.value_counts().mean())
        # print(1 + taxa.value_counts().std())

        # num_clusters_modifier = 1 - abs(len(self.clusters) - genes_per_taxa) / genes_per_taxa

        # modifier = abs(taxa.value_counts().mean() ** (1 + taxa.value_counts().std()) - len(self.clusters))
        # sys.exit("%s" % (modifier))
        # return (score + num_clusters_modifier) / 2
        return self.srs(score)

    # The following are a group of possible scoring schemes
    def gpt(self, score, taxa):  # groups per taxa
        genes_per_taxa = self.size / len(taxa.value_counts())
        num_clusters_modifier = abs(genes_per_taxa - len(self.clusters))
        score = round(score / len(self.clusters) - num_clusters_modifier ** 1.2, 2)
        return score

    def srs(self, score):  # Square root sequences
        # print(len(self.clusters))
        modifier = abs(self.size ** (1 / 2) - len(self.clusters)) ** 1.2
        score = round(score / len(self.clusters) - modifier, 2)
        return score

    def write(self, outfile):
        with open(outfile, "w") as _ofile:
            for cluster in self.clusters:
                _ofile.write("%s\n" % "\t".join(cluster.cluster))
        return


class Cluster:
    def __init__(self, cluster):
        cluster.sort()
        self.cluster = cluster
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

        # perfect_score = len(taxa) ** 2
        # score += singles ** 2
        # score /= perfect_score
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
    external_tmp_dir, min_score = params
    mcl_tmp_dir = MyFuncs.TempDir()

    _output = Popen("mcl %s/input.csv --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" %
                    (external_tmp_dir, gq, I, mcl_tmp_dir.path), shell=True, stderr=PIPE).communicate()

    _output = _output[1].decode()

    if re.search("\[mclvInflate\] warning", _output) and min_score:
        return min_score

    clusters = Clusters("%s/output.groups" % mcl_tmp_dir.path)
    score = clusters.score_all_clusters()

    with lock:
        with open("%s/max.txt" % external_tmp_dir, "r") as _ifile:
            current_max = float(_ifile.read())

    if score > current_max:
        with lock:
            clusters.write("%s/output.groups" % external_tmp_dir)
            with open("%s/max.txt" % external_tmp_dir, "w") as _ofile:
                _ofile.write(str(score))

    return score


def homolog_caller(cluster, all_by_all, cluster_list, rank, save=False, steps=1000):
    temp_dir = MyFuncs.TempDir()

    all_by_all.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    inflation_var = mcmcmc.Variable("I", 1.4, 20)
    gq_var = mcmcmc.Variable("gq", 0.05, 0.95)

    try:
        with open("%s/max.txt" % temp_dir.path, "w") as _ofile:
            _ofile.write("-1000000000")

        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1,
                                       params=["%s" % temp_dir.path, False], quiet=True,
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

    mcmcmc_factory.reset_params(["%s" % temp_dir.path, worst_score])

    mcmcmc_factory.run()

    mcmcmc_output = pd.read_csv("%s/mcmcmc_out.csv" % temp_dir.path, "\t")

    best_score = max(mcmcmc_output["result"])
    if best_score <= cluster.score():
        cluster_list.append(cluster)
        if save:
            temp_dir.save("%s/group_%s" % (save, rank))
        return cluster_list

    mcl_clusters = Clusters("%s/output.groups" % temp_dir.path)
    _counter = 1
    for _clust in mcl_clusters.clusters:
        next_rank = "%s_%s" % (rank, _counter)
        _clust.name = next_rank
        _counter += 1
        if len(_clust.cluster) in [1, 2]:
            cluster_list.append(_clust)
            continue

        group_all_by_all = split_all_by_all(all_by_all, _clust.cluster)["removed"]

        _min = group_all_by_all[:][2].min()
        data_range = group_all_by_all[:][2].max() - _min

        re_norm = (group_all_by_all[:][2] - _min) / data_range

        group_all_by_all[:][2] = re_norm

        # Recursion...
        cluster_list = homolog_caller(_clust, group_all_by_all, cluster_list, next_rank, save, steps=steps)

    if save:
        temp_dir.save("%s/group_%s" % (save, rank))

    return cluster_list


def merge_singles(clusters, scores):
    small_clusters = []
    large_clusters = []
    large_group_names = []
    small_group_names = []
    for cluster in clusters:
        if len(cluster.cluster) > 2:
            large_group_names.append(cluster.name)
            large_clusters.append(cluster)
        else:
            small_group_names.append(cluster.name)
            small_clusters.append(cluster)

    # Convert the large_clusters list to a dict using group name as key
    large_clusters = {x: large_clusters[j] for j, x in enumerate(large_group_names)}

    small_to_large_dict = {}
    for sclust in small_clusters:
        small_to_large_dict[sclust.name] = {ind: [] for ind, value in large_clusters.items()}
        for sgene in sclust.cluster:
            for key, lclust in large_clusters.items():
                for lgene in lclust.cluster:
                    score = scores.loc[:][scores[0] == sgene]
                    score = score.loc[:][score[1] == lgene]

                    if score.empty:
                        score = scores.loc[:][scores[0] == lgene]
                        score = score.loc[:][score[1] == sgene]

                    score = float(score[2])
                    small_to_large_dict[sclust.name][lclust.name].append(score)

    small_clusters = {x: small_clusters[j] for j, x in enumerate(small_group_names)}
    for small_group_id, l_clusts in small_to_large_dict.items():
        # Convert data into list of numpy arrays that sm.stats can read, also get average scores for each cluster
        data = [np.array(x) for x in l_clusts.values()]
        averages = pd.Series()
        for j, group in enumerate(data):
            key = list(l_clusts.keys())[j]
            averages = averages.append(pd.Series(np.mean(group), index=[key]))
            data[j] = pd.DataFrame(group, columns=['observations'])
            data[j]['grouplabel'] = key

        max_ave = averages.argmax()

        df = pd.DataFrame()
        for group in data:
            df = df.append(group)

        result = sm.stats.multicomp.pairwise_tukeyhsd(df.observations, df.grouplabel)

        for line in str(result).split("\n")[4:-1]:
            line = re.sub("^ *", "", line.strip())
            line = re.sub(" +", "\t", line)
            line = line.split("\t")
            if max_ave in line:
                if line[5] == 'True':
                    continue
                else:
                    # Insufficient support to group the gene with max_ave group
                    break
        else:
            # The gene can be grouped with the max_ave group (write the code)
            # final_clusters[max_ave].cluster.append(small_group_id)
            large_clusters[max_ave].cluster += small_clusters[small_group_id].cluster
            del small_clusters[small_group_id]

    clusters = [value for ind, value in large_clusters.items()]
    clusters += [value for ind, value in small_clusters.items()]
    return clusters


def jackknife(orig_clusters, all_by_all, steps=1000, level=0.632):
    total_pop = [x.cluster for x in orig_clusters]  # List of lists
    total_pop = [x for _clust in total_pop for x in _clust]  # Flatten to a 1D list

    sample_size = floor(level * len(total_pop))

    # Add a new 'support' attribute to each Cluster object to be filled in with each step
    for _clust in orig_clusters:
        _clust.support = 0.

    for _ in range(steps):
        jn_sample = Cluster(sample(total_pop, sample_size))
        samp_all_by_all = split_all_by_all(all_by_all, jn_sample.cluster)["remaining"]
        sample_clusters = []
        sample_clusters = homolog_caller(jn_sample, samp_all_by_all, sample_clusters, 0, False, steps=1000)
        sample_clusters = merge_singles(sample_clusters, samp_all_by_all)

        for orig_clust in orig_clusters:
            for _clust in sample_clusters:
                intersect = set(_clust.cluster).intersection(orig_clust.cluster)
                if len(intersect) == 0:
                    continue

                elif len(intersect) != len(_clust.cluster):
                    break

                else:
                    orig_clust.support += 1

        ####
        _output = ""
        for _clust in sample_clusters:
            _output += "group_%s\t%s\t" % (_clust.name, _clust.score())
            for _seq_id in _clust.cluster:
                _output += "%s\t" % _seq_id
            _output = "%s\n" % _output.strip()
        print(_output)
        ######

    for _clust in orig_clusters:
        _clust.support /= steps

    return orig_clusters

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="homolog_caller", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("all_by_all_file", help="Location of all-by-all scores file", action="store")
    parser.add_argument("output_file", help="Where should groups be written to? Default to stdout.", action="store", nargs="?")
    parser.add_argument("-jk", "--jackknife", help="Find support for previous run.", metavar="<Groups file>")
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

    if in_args.jackknife:
        print("Preparing Jackknife support values")
        with open(in_args.jackknife, "r") as ifile:
            groups = ifile.read()

        groups = groups.strip().split("\n")
        final_clusters = [group.strip().split("\t") for group in groups]
        for i, clust in enumerate(final_clusters):
            final_clusters[i] = Cluster(clust[2:])
            final_clusters[i].name = clust[0]

        final_clusters = jackknife(final_clusters, scores_data, 100)
        for clust in final_clusters:
            print("%s: %s" % (clust.name, clust.support))

    else:
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

        # Try to fold singletons and doublets back into groups.
        final_clusters = merge_singles(final_clusters, scores_data)

        output = ""
        while len(final_clusters) > 0:
            _max = (0, 0)
            for ind, clust in enumerate(final_clusters):
                if len(clust.cluster) > _max[1]:
                    _max = (ind, len(clust.cluster))

            ind, _max = _max[0], final_clusters[_max[0]]
            del final_clusters[ind]
            output += "group_%s\t%s\t" % (_max.name, _max.score())
            for seq_id in _max.cluster:
                output += "%s\t" % seq_id
            output = "%s\n" % output.strip()

        if not in_args.output_file:
            print(output)

        else:
            with open(in_args.output_file, "w") as ofile:
                ofile.write(output)
