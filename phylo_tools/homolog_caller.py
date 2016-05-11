#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 27 2015 

"""
Convert an all-by-all scores matrix into groups. MCMCMC-driven MCL is the basic work horse for clustering, and then
the groups are refined further to include singletons and doublets and/or be broken up into cliques, where appropriate.
"""

# Std library
import sys
import os
import re
import shutil
from copy import copy
from subprocess import Popen, PIPE
from multiprocessing import Lock
from random import random

# 3rd party
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats
from Bio.SubsMat import SeqMat, MatrixInfo

# My packages
import mcmcmc  # Note: This is in ../utilities and sym-linked to python3.5/site-packages
import MyFuncs
import SeqBuddy as Sb
import AlignBuddy as Alb


class Clusters:  # The cluster groups should not include anything more than the groups (no names)
    def __init__(self, path, group_split="\n", gene_split="\t", taxa_split="-", global_taxa_count=None):
        with open(path, "r") as _ifile:
            self.input = _ifile.read()

        clusters = self.input.strip().split(group_split)
        self.clusters = [Cluster([y for y in x.strip().split(gene_split)], global_taxa_count=global_taxa_count)
                         for x in clusters]
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
            score += cluster.score() * (cluster.len / self.size)
            temp_taxa = [x.split(self.taxa_split)[0] for x in cluster.cluster]
            temp_taxa = pd.Series(temp_taxa)
            taxa = pd.concat([taxa, temp_taxa])

        return self.srs(score)

    def write(self, outfile):
        with open(outfile, "w") as _ofile:
            for cluster in self.clusters:
                _ofile.write("%s\n" % "\t".join(cluster.cluster))
        return

    # The following are possible score modifier schemes to account for group size
    def gpt(self, score, taxa):  # groups per taxa
        genes_per_taxa = self.size / len(taxa.value_counts())
        num_clusters_modifier = abs(genes_per_taxa - len(self.clusters))
        score = round(score / len(self.clusters) - num_clusters_modifier ** 1.2, 2)
        return score

    def srs(self, score):  # Square root sequences
        # print(len(self.clusters))
        modifier = abs(self.size ** 0.5 - len(self.clusters)) ** 1.2
        score = round(score / len(self.clusters) - modifier, 2)
        return score


class Cluster:
    def __init__(self, cluster, taxa_split="-", global_taxa_count=None, _name=""):
        cluster.sort()
        self.cluster = cluster
        self.len = len(self.cluster)
        self.name = _name
        self.taxa_split = taxa_split
        self.global_taxa_count = global_taxa_count

    def score(self):
        taxa = [x.split(self.taxa_split)[0] for x in self.cluster]
        taxa = pd.Series(taxa)

        if len(taxa) == 1:
            return 0

        score = 0
        try:
            for taxon, num in taxa.value_counts().iteritems():
                if num == 1:
                    score += self.global_taxa_count[taxon] ** 0.5

                else:
                    score -= (1 / self.global_taxa_count[taxon] ** 0.5) * num * 2

        except TypeError:  # This happens if global_taxa_count is not set
            singles = 0
            for taxon, num in taxa.value_counts().iteritems():
                if num == 1:
                    singles += 1

                else:
                    score -= num

            score += singles ** 2

        return score

    def compare(self, query):
        matches = set(self.cluster).intersection(query.cluster)
        weighted_match = (len(matches) * 2.) / (self.len + query.len)
        print("name: %s, matches: %s, weighted_match: %s" % (self.name, len(matches), weighted_match))
        return weighted_match

    def __str__(self):
        return str(self.cluster)


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
    inflation, gq = args
    external_tmp_dir, min_score, all_by_all, global_taxa_count = params
    mcl_tmp_dir = MyFuncs.TempDir()

    _output = Popen("mcl %s/input.csv --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" %
                    (external_tmp_dir, gq, inflation, mcl_tmp_dir.path), shell=True, stderr=PIPE).communicate()

    _output = _output[1].decode()

    if re.search("\[mclvInflate\] warning", _output) and min_score:
        return min_score

    clusters = Clusters("%s/output.groups" % mcl_tmp_dir.path, global_taxa_count=global_taxa_count)
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


def clique_checker(cluster, df_all_by_all):
    def get_best_hit(gene_name):  # ToDo: When pulling in best hits, check for multiple equal best hits
        _best_hit = df_all_by_all[(df_all_by_all[0] == gene_name) | (df_all_by_all[1] == gene_name)]
        _best_hit.columns = ["subj", "query", "score"]
        _best_hit = _best_hit.loc[_best_hit["score"] == max(_best_hit["score"])].values[0]
        return _best_hit

    valve = MyFuncs.SafetyValve(50)
    genes = pd.DataFrame([x.split(cluster.taxa_split) for x in cluster.cluster])
    genes.columns = ["taxa", "gene"]

    in_dict = {}
    cliques = []

    # pull out all genes with replicate taxa and get best hits
    for taxa, num in genes['taxa'].value_counts().iteritems():
        if num > 1:
            taxa = genes.loc[genes["taxa"] == taxa]
            for j, rec in taxa.iterrows():
                _gene = "%s-%s" % (rec["taxa"], rec["gene"])
                best_hit = get_best_hit(_gene)
                # Don't know what order the "subj" and "query" genes are in, so check
                # in_dict[_gene] = (<the best hit gene>, <best hit sim score>)
                in_dict[_gene] = (best_hit[0], best_hit[2]) if best_hit[1] == _gene else (best_hit[1], best_hit[2])

    # Iterate through in_dict and pull in all genes to fill out all best-hit sub-graphs from duplicates
    in_dict_length = len(in_dict)
    while True:
        copy_in_dict = copy(in_dict)
        valve.step()
        for indx, value in copy_in_dict.items():
            _gene = value[0]
            if _gene not in in_dict:
                best_hit = get_best_hit(_gene)
                in_dict[_gene] = (best_hit[0], best_hit[2]) if best_hit[1] == _gene else (best_hit[1], best_hit[2])

        if in_dict_length == len(in_dict):
            break
        else:
            in_dict_length = len(in_dict)

    # Separate all sub-graphs (i.e., cliques) into lists
    while len(in_dict):
        valve.step()
        copy_in_dict = copy(in_dict)
        for indx, value in copy_in_dict.items():
            in_clique = False

            for j, clique in enumerate(cliques):
                if indx in clique or value[0] in clique:
                    if type(in_clique) != int:
                        in_clique = j
                        cliques[j] = list({indx, value[0]}.union(set(cliques[j])))
                        del in_dict[indx]

                    else:
                        cliques[in_clique] = list(set(cliques[in_clique]).union(set(cliques[j])))
                        del cliques[j]
                        break

            if type(in_clique) != int:
                del in_dict[indx]
                cliques.append([indx, value[0]])

    # Strip out any 'cliques' that contain less than 3 genes
    while True:
        valve.step()
        for j, clique in enumerate(cliques):
            if len(clique) < 3:
                del cliques[j]
                break
        break

    # Get the similarity scores for within cliques and between cliques-remaining sequences, then generate kernel-density
    # functions for both. The overlap between the two functions is used to determine whether they should be separated
    final_cliques = []
    for clique in cliques:
        total_scores = pd.DataFrame()
        if len(clique) < 3:
            continue

        for _gene in clique:
            scores = df_all_by_all[df_all_by_all[0] == _gene]
            scores = scores[scores[1].isin(cluster.cluster)]
            tmp = df_all_by_all[df_all_by_all[1] == _gene]
            tmp = tmp[tmp[0].isin(cluster.cluster)]
            scores = pd.concat([scores, tmp])
            total_scores = total_scores.append(scores)

        clique_scores = total_scores[(total_scores[0].isin(clique)) & (total_scores[1].isin(clique))]
        total_scores = total_scores.drop(clique_scores.index.values)

        # if a clique is found that pulls in every single gene, skip
        if not len(total_scores):
            continue

        # if all sim scores in a group are identical, we can't get a KDE. Fix by perturbing the scores a little.
        if clique_scores[2].std() == 0:
            for indx, score in clique_scores[2].iteritems():
                min_max = [score - (score * 0.01), score + (score * 0.01)]
                change = random() * (min_max[1] - min_max[0])
                clique_scores.loc[indx, 2] = min_max[0] + change
        if total_scores[2].std() == 0:
            for indx, score in total_scores[2].iteritems():
                min_max = [score - (score * 0.01), score + (score * 0.01)]
                change = random() * (min_max[1] - min_max[0])
                total_scores.loc[indx, 2] = min_max[0] + change

        total_kde = scipy.stats.gaussian_kde(total_scores[2], bw_method='silverman')
        clique_kde = scipy.stats.gaussian_kde(clique_scores[2], bw_method='silverman')
        clique_resample = clique_kde.resample(10000)
        clique95 = [scipy.stats.scoreatpercentile(clique_resample, 2.5),
                    scipy.stats.scoreatpercentile(clique_resample, 97.5)]

        integrated = total_kde.integrate_box_1d(clique95[0], clique95[1])
        if integrated < 0.05:
            final_cliques.append(clique)

    if final_cliques:
        for clique in final_cliques:
            cluster.cluster = [x for x in cluster.cluster if x not in clique]

        final_cliques = [Cluster(x) for x in final_cliques]
        for j, clique in enumerate(final_cliques):
            clique.name = "%s_c%s" % (cluster.name, j + 1)

    final_cliques.append(cluster)
    return final_cliques


def homolog_caller(cluster, local_all_by_all, cluster_list, rank, seqbuddy, global_all_by_all=None, steps=1000,
                   global_taxa_count=None, quiet=True, _clique_check=True, recursion=True, ):

    temp_dir = MyFuncs.TempDir()

    local_all_by_all.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    inflation_var = mcmcmc.Variable("I", 1.1, 20)
    gq_var = mcmcmc.Variable("gq", min(local_all_by_all[2]), max(local_all_by_all[2]))

    try:
        with open("%s/max.txt" % temp_dir.path, "w") as _ofile:
            _ofile.write("-1000000000")

        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1,
                                       params=["%s" % temp_dir.path, False, False, global_taxa_count], quiet=quiet,
                                       outfile="%s/mcmcmc_out.csv" % temp_dir.path)

    except RuntimeError:  # Happens when mcmcmc fails to find different initial chain parameters
        cluster_list.append(cluster)
        temp_dir.save("%s/mcmcmc/group_%s" % (in_args.outdir, rank))
        return cluster_list

    # Set a 'worst score' that is reasonable for the data set
    worst_score = 10000000  # arbitrarily large number to start
    for chain in mcmcmc_factory.chains:
        worst_score = chain.raw_min if chain.raw_min < worst_score else worst_score

    mcmcmc_factory.reset_params(["%s" % temp_dir.path, worst_score, global_all_by_all, global_taxa_count])
    mcmcmc_factory.run()
    mcmcmc_output = pd.read_csv("%s/mcmcmc_out.csv" % temp_dir.path, "\t")

    best_score = max(mcmcmc_output["result"])
    if best_score <= cluster.score():
        if _clique_check:
            cluster_list += clique_checker(cluster, local_all_by_all)
        else:
            cluster_list.append(cluster)
        temp_dir.save("%s/mcmcmc/group_%s" % (in_args.outdir, rank))
        return cluster_list

    mcl_clusters = Clusters("%s/output.groups" % temp_dir.path, global_taxa_count=global_taxa_count)
    _counter = 1
    for _clust in mcl_clusters.clusters:
        next_rank = "%s_%s" % (rank, _counter)
        _clust.name = next_rank
        _counter += 1
        if len(_clust.cluster) in [1, 2]:
            cluster_list.append(_clust)
            continue

        seqbuddy_copy = Sb.make_copy(seqbuddy)
        seqbuddy_copy = Sb.pull_recs(seqbuddy_copy, ["^%s$" % rec_id for rec_id in _clust.cluster])
        _scores_data = create_all_by_all_scores(seqbuddy_copy, group=_clust.name)

        sub_cluster = pd.concat([_scores_data[0], _scores_data[1]])
        sub_cluster = sub_cluster.value_counts()
        sub_cluster = Cluster([i for i in sub_cluster.index], _name=next_rank)

        # Recursion...
        if recursion:
            cluster_list = homolog_caller(sub_cluster, _scores_data, cluster_list, next_rank, seqbuddy=seqbuddy_copy,
                                          steps=steps, global_all_by_all=global_all_by_all,
                                          global_taxa_count=global_taxa_count, quiet=quiet, _clique_check=_clique_check)
        else:
            cluster_list.append(_clust)

    temp_dir.save("%s/mcmcmc/group_%s" % (in_args.outdir, rank))
    return cluster_list


def merge_singles(clusters, scores):
    small_clusters = []
    small_group_names = []
    large_clusters = []
    large_group_names = []
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
        small_to_large_dict[sclust.name] = {_ind: [] for _ind, value in large_clusters.items()}
        for sgene in sclust.cluster:
            for key, lclust in large_clusters.items():
                for lgene in lclust.cluster:
                    score = scores.loc[:][scores[0] == sgene]
                    score = score.loc[:][score[1] == lgene]

                    if score.empty:
                        score = scores.loc[:][scores[0] == lgene]
                        score = score.loc[:][score[1] == sgene]

                    if score.empty:
                        score = 0.
                    else:
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
            # The gene can be grouped with the max_ave group
            large_clusters[max_ave].cluster += small_clusters[small_group_id].cluster
            del small_clusters[small_group_id]

    clusters = [value for _ind, value in large_clusters.items()]
    clusters += [value for _ind, value in small_clusters.items()]
    return clusters
# NOTE: There used to be a support function. Check the GitHub history if there's desire to bring parts of it back


def score_sequences(pair, args):
    id1, id2 = pair
    alignbuddy, outfile = args
    id_regex = "^%s$|^%s$" % (id1, id2)
    alb_copy = Alb.make_copy(alignbuddy)
    Alb.pull_records(alb_copy, id_regex)
    _score = 0
    seq1, seq2 = alb_copy.records()
    prev_aa1 = "-"
    prev_aa2 = "-"

    for aa_pos in range(alb_copy.lengths()[0]):
        aa1 = seq1.seq[aa_pos]
        aa2 = seq2.seq[aa_pos]

        if aa1 == "-" or aa2 == "-":
            if prev_aa1 == "-" or prev_aa2 == "-":
                _score += gap_extend
            else:
                _score += gap_open
        else:
            _score += BLOSUM62[aa1, aa2]
        prev_aa1 = str(aa1)
        prev_aa2 = str(aa2)

    with lock:
        with open(outfile, "a") as _ofile:
            _ofile.write("%s\t%s\t%s\n" % (id1, id2, _score))
    return


def create_all_by_all_scores(seqs, group, progress=False):
    alignment = Alb.generate_msa(Sb.make_copy(seqs), tool="mafft", params="--globalpair", quiet=True)
    alignment = Alb.trimal(alignment, "gappyout")
    alignment.write("%s/alignments/group_%s.aln" % (in_args.outdir, group))
    ids1 = [rec.id for rec in alignment.records_iter()]
    ids2 = [rec.id for rec in alignment.records_iter()]
    all_by_all = []
    for rec1 in ids1:
        del ids2[ids2.index(rec1)]
        for rec2 in ids2:
            all_by_all.append((rec1, rec2))
    outfile = "%s/sim_scores/group_%s.csv" % (in_args.outdir, group)
    MyFuncs.run_multicore_function(all_by_all, score_sequences, [alignment, outfile])
    _output = pd.read_csv(outfile, sep="\t", header=None)
    return _output

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="homolog_caller", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("sequences", help="Location of sequence file", action="store")
    parser.add_argument("outdir", action="store", default="%s/rd-mcd" % os.getcwd(),
                        help="Where should results be written?")
    parser.add_argument("-sz", "--sample_size", type=float, default=0.632,
                        help="Proportion of total population to use in each jackknife replicate")
    parser.add_argument("-mcs", "--mcmcmc_steps", default=1000, type=int,
                        help="Specify how deeply to sample MCL parameters")
    parser.add_argument("-sep", "--separator", action="store", default="\t",
                        help="If the all-by-all file is not tab delimited, specify the character")
    parser.add_argument("-sr", "--supress_recursion", action="store_true",
                        help="Stop after a single round of MCL. For testing.")
    parser.add_argument("-scc", "--supress_clique_check", action="store_true",
                        help="Do not check for or break up cliques. For testing.")
    parser.add_argument("-ssf", "--supress_singlet_folding", action="store_true",
                        help="Do not check for or merge singlets. For testing.")
    parser.add_argument("-op", "--open_penalty", help="Penalty for opening a gap in pairwise alignment scoring",
                        type=float, default=-10)
    parser.add_argument("-ep", "--extend_penalty", help="Penalty for extending a gap in pairwise alignment scoring",
                        type=float, default=-1)
    parser.add_argument("-nt", "--no_msa_trim", action="store_true",
                        help="Don't apply the gappyout algorithm to MSAs before scoring")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Overwrite previous run")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress all output during run (only final output is returned)")

    in_args = parser.parse_args()

    if os.path.exists(in_args.outdir):
        check = MyFuncs.ask("Output directory already exists, overwrite it [y]/n?") if not in_args.force else True
        if check:
            shutil.rmtree(in_args.outdir)
        else:
            print("Program aborted. Output directory required.")
            sys.exit()

    clique_check = True if not in_args.supress_clique_check else False
    recursion_check = True if not in_args.supress_recursion else False
    best = None
    best_clusters = None
    lock = Lock()
    printer = MyFuncs.DynamicPrint()

    sequences = Sb.SeqBuddy(in_args.sequences)
    PHAT = make_full_mat(SeqMat(MatrixInfo.phat75_73))
    BLOSUM62 = make_full_mat(SeqMat(MatrixInfo.blosum62))
    BLOSUM45 = make_full_mat(SeqMat(MatrixInfo.blosum45))

    gap_open = in_args.open_penalty
    gap_extend = in_args.extend_penalty

    os.makedirs(in_args.outdir)
    os.makedirs("%s/alignments" % in_args.outdir)
    os.makedirs("%s/mcmcmc" % in_args.outdir)
    os.makedirs("%s/sim_scores" % in_args.outdir)

    print("Generating initial all-by-all...")
    scores_data = create_all_by_all_scores(sequences, group="0", progress=True)

    master_cluster = pd.concat([scores_data[0], scores_data[1]])
    master_cluster = master_cluster.value_counts()
    master_cluster = Cluster([i for i in master_cluster.index], _name="0")

    taxa_count = [x.split("-")[0] for x in master_cluster.cluster]
    taxa_count = pd.Series(taxa_count)
    taxa_count = taxa_count.value_counts()

    print("Creating clusters")
    final_clusters = []
    final_clusters = homolog_caller(master_cluster, scores_data, final_clusters, rank=0, seqbuddy=sequences,
                                    global_all_by_all=scores_data,
                                    steps=in_args.mcmcmc_steps, global_taxa_count=taxa_count, quiet=False,
                                    _clique_check=clique_check, recursion=recursion_check)

    # Try to fold singletons and doublets back into groups.
    if not in_args.supress_singlet_folding:
        final_clusters = merge_singles(final_clusters, scores_data)

    # Format the clusters and output to stdout or file
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

    with open("%s/clusters.txt" % in_args.outdir, "w") as ofile:
        ofile.write(output)
