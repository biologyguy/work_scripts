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
from math import log

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


class ClusterNew(object):
    def __init__(self, cluster, sim_scores, _parent=None, clique=False):
        """
        - Note that reciprocal best hits between paralogs are collapsed when instantiating group_0, so
          no problem strongly penalizing all paralogs in the scoring algorithm

        :param cluster: Sequence IDs
        :type cluster: list
        :param sim_scores: All-by-all similarity matrix for everything in the cluster (and parental clusters)
        :type sim_scores: pandas.DataFrame
        :param _parent: Parental cluster
        :type _parent: Cluster
        """
        self.taxa = {}
        self.sim_scores = sim_scores
        self.parent = _parent
        self._subgroup_counter = 0
        self.cliques = None
        self._score = None
        self.collapsed_genes = {}  # If paralogs are reciprocal best hits, collapse them

        if clique and not _parent:
            raise AttributeError("A clique cannot be declared without including its parental cluster.")

        if _parent:
            if clique:
                self.name = "%s_c0" % _parent._name if not _parent.cliques else \
                    "%s_c%s" % (_parent._name, len(_parent.cliques) + 1)
            else:
                self._name = "%s_%s" % (_parent._name, _parent._subgroup_counter)
                _parent._subgroup_counter += 1
            for indx, genes in _parent.collapsed_genes.items():
                if indx in cluster:
                    self.collapsed_genes[indx] = genes
            self.cluster = cluster
            for gene in cluster:
                taxa = gene.split("-")[0]
                self.taxa.setdefault(taxa, [])
                self.taxa[taxa].append(gene)
        else:
            self._name = "group_0"
            collapse_check_list = []
            collapsed_cluster = []
            for gene in cluster:
                if gene in collapse_check_list:
                    continue

                taxa = gene.split("-")[0]
                self.taxa.setdefault(taxa, [])
                self.taxa[taxa].append(gene)
                collapsed_cluster.append(gene)
                valve = MyFuncs.SafetyValve()  # ToDo: remove this when sure it is unnecessary
                breakout = False
                while not breakout:  # This is the primary paralog collapsing logic. Only do it once for parental group
                    valve.step()
                    breakout = True
                    for edge in self._get_best_hits(gene).itertuples():
                        other_seq_id = edge.seq1 if edge.seq1 != gene else edge.seq2
                        if other_seq_id.split("-")[0] == taxa:
                            other_seq_best_hits = self._get_best_hits(other_seq_id)
                            if gene in other_seq_best_hits.seq1.values or gene in other_seq_best_hits.seq2.values:
                                self.collapsed_genes.setdefault(gene, [])
                                self.collapsed_genes[gene].append(other_seq_id)
                                collapse_check_list.append(other_seq_id)
                                # Strip collapsed paralogs from the all-by-all graph
                                self.sim_scores = self.sim_scores[(self.sim_scores.seq1 != other_seq_id) &
                                                                  (self.sim_scores.seq2 != other_seq_id)]
                                if other_seq_id in collapsed_cluster:
                                    del collapsed_cluster[collapsed_cluster.index(other_seq_id)]
                                    del self.taxa[taxa][self.taxa[taxa].index(other_seq_id)]
                                    if other_seq_id in self.collapsed_genes:
                                        collapse_check_list += self.collapsed_genes[other_seq_id]
                                        del self.collapsed_genes[other_seq_id]
                                breakout = False
                                break
            self.cluster = collapsed_cluster

    def _get_best_hits(self, gene):
        best_hits = self.sim_scores[(self.sim_scores.seq1 == gene) | (self.sim_scores.seq2 == gene)]
        try:
            best_hits = best_hits.loc[best_hits.score == max(best_hits.score)].values
        except ValueError:
            print(gene)
            sys.exit()
        best_hits = pd.DataFrame(best_hits, columns=["seq1", "seq2", "score"])
        return best_hits

    def compare(self, query):
        matches = set(self.cluster).intersection(query.cluster)
        weighted_match = (len(matches) * 2.) / (len(self) + query.len)
        print("name: %s, matches: %s, weighted_match: %s" % (self.name, len(matches), weighted_match))
        return weighted_match

    def sub_cluster(self, id_list):
        sim_scores = pd.DataFrame()
        for gene in id_list:
            if gene not in self.cluster:
                raise NameError("Gene id '%s' not present in parental cluster" % gene)
            seq1_scores = self.sim_scores[self.sim_scores.seq1 == gene]
            seq1_scores = seq1_scores[seq1_scores.seq2.isin(id_list)]
            seq2_scores = self.sim_scores[self.sim_scores.seq2 == gene]
            seq2_scores = seq2_scores[seq2_scores.seq1.isin(id_list)]
            scores = pd.concat([seq1_scores, seq2_scores])
            sim_scores = sim_scores.append(scores)
        sim_scores = sim_scores.drop_duplicates()
        subcluster = Cluster(id_list, sim_scores, self)
        return subcluster

    def name(self):
        return self._name

    def _recursive_best_hits(self, gene, global_best_hits, tested_ids):
        best_hits = self._get_best_hits(gene)
        global_best_hits = global_best_hits.append(best_hits, ignore_index=True)
        for _edge in best_hits.itertuples():
            if _edge.seq1 not in tested_ids:
                tested_ids.append(_edge.seq1)
                global_best_hits = self._recursive_best_hits(_edge.seq1, global_best_hits, tested_ids)
            if _edge.seq2 not in tested_ids:
                tested_ids.append(_edge.seq2)
                global_best_hits = self._recursive_best_hits(_edge.seq2, global_best_hits, tested_ids)
        return global_best_hits

    def clique_checker(self):
        best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
        # pull out all genes with replicate taxa and get best hits
        for taxa, genes in self.taxa.items():
            if len(genes) > 1:
                for gene in genes:
                    best_hits = self._recursive_best_hits(gene, best_hits, [gene])

        cliques = []
        for edge in best_hits.itertuples():
            match_indicies = []
            for indx, clique in enumerate(cliques):
                if edge.seq1 in clique and indx not in match_indicies:
                    match_indicies.append(indx)
                if edge.seq2 in clique and indx not in match_indicies:
                    match_indicies.append(indx)
                if len(match_indicies) == 2:
                    break

            if not match_indicies:
                cliques.append([edge.seq1, edge.seq2])
            elif len(match_indicies) == 1:
                new_clique = set(cliques[match_indicies[0]] + [edge.seq1, edge.seq2])
                cliques[match_indicies[0]] = list(new_clique)
            else:
                match_indicies.sort()
                new_clique = set(cliques[match_indicies[0]] + cliques[match_indicies[1]])
                cliques[match_indicies[0]] = list(new_clique)
                del cliques[match_indicies[1]]

        # Strip out any 'cliques' that contain less than 3 genes
        cliques = [clique for clique in cliques if len(clique) >= 3]

        # Get the similarity scores for within cliques and between cliques-remaining sequences, then generate kernel-density
        # functions for both. The overlap between the two functions is used to determine whether they should be separated

        def perturb(_scores):
            _valve = MyFuncs.SafetyValve(global_reps=10)
            while _scores.score.std() == 0:
                _valve.step("Failed to perturb:\n%s" % _scores)
                for _indx, _score in _scores.score.iteritems():
                    _scores.set_value(_indx, "score", random.gauss(_score, (_score * 0.0000001)))
            return _scores

        self.cliques = []
        for clique in cliques:
            clique_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique)) & (self.sim_scores.seq2.isin(clique))]
            total_scores = self.sim_scores.drop(clique_scores.index.values)

            # if a clique is found that pulls in every single gene, skip
            if not len(total_scores):
                continue

            # if all sim scores in a group are identical, we can't get a KDE. Fix by perturbing the scores a little.
            clique_scores = perturb(clique_scores)
            total_scores = perturb(total_scores)

            total_kde = scipy.stats.gaussian_kde(total_scores.score, bw_method='silverman')
            clique_kde = scipy.stats.gaussian_kde(clique_scores.score, bw_method='silverman')
            clique_resample = clique_kde.resample(10000)
            clique95 = [scipy.stats.scoreatpercentile(clique_resample, 2.5),
                        scipy.stats.scoreatpercentile(clique_resample, 97.5)]

            integrated = total_kde.integrate_box_1d(clique95[0], clique95[1])
            if integrated < 0.05:
                clique = Cluster(clique, sim_scores=clique_scores, _parent=self, clique=True)
                self.cliques.append(clique)
        self.cliques = [None] if not self.cliques else self.cliques
        return

    def score(self):
        if self._score:
            return self._score
        # Don't ignore the possibility of cliques, which will alter the score.
        if not self.cliques:
            self.clique_checker()
        if self.cliques[0]:
            clique_list = [i for j in self.cliques for i in j.cluster]
            decliqued_cluster = []
            for gene in self.cluster:
                if gene not in clique_list:
                    decliqued_cluster.append(gene)
            score = self.raw_score(decliqued_cluster)
            score += sum([self.raw_score(x.cluster) for x in self.cliques])
        else:
            score = self.raw_score(self.cluster)
        self._score = score
        return score

    def raw_score(self, id_list):
        if len(id_list) == 1:
            return 0

        taxa = {}
        for gene in id_list:
            taxon = gene.split("-")[0]
            taxa.setdefault(taxon, [])
            taxa[taxon].append(gene)

        # An entire cluster should never consist of a single taxa because I've stripped out reciprocal best hits paralogs
        if len(taxa) == 1:
            raise ReferenceError("Only a single taxa found in cluster...")

        unique_scores = 1
        paralog_scores = -1
        for taxon, genes in taxa.items():
            if len(genes) == 1:
                unique_scores *= 2 * (1 + (len(self.taxa[taxon]) / len(self)))
            else:
                paralog_scores *= len(genes) * (len(genes) / len(self.taxa[taxon]))
        return unique_scores + paralog_scores

    def __len__(self):
        return len(self.cluster)

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


def _psi_pred(seq_obj):
    if os.path.isfile("%s/psi_pred/%s.ss2" % (in_args.outdir, seq_obj.id)):
        return
    temp_dir = MyFuncs.TempDir()
    pwd = os.getcwd()
    os.chdir(temp_dir.path)
    with open("sequence.fa", "w") as _ofile:
        _ofile.write(seq_obj.format("fasta"))

    Popen("runpsipred sequence.fa > /dev/null 2>&1", shell=True).wait()
    os.chdir(pwd)
    shutil.move("%s/sequence.ss2" % temp_dir.path, "%s/psi_pred/%s.ss2" % (in_args.outdir, seq_obj.id))
    return


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


def score_sequences(_pair, args):
    # Calculate the best possible scores, and divide by the observed scores
    id1, id2 = _pair
    alignbuddy, psi_pred_files, outfile = args
    id_regex = "^%s$|^%s$" % (id1, id2)
    alb_copy = Alb.make_copy(alignbuddy)
    Alb.pull_records(alb_copy, id_regex)
    observed_score = 0
    seq1_best = 0
    seq2_best = 0
    seq1, seq2 = alb_copy.records()
    prev_aa1 = "-"
    prev_aa2 = "-"

    for aa_pos in range(alb_copy.lengths()[0]):
        aa1 = seq1.seq[aa_pos]
        aa2 = seq2.seq[aa_pos]

        if aa1 != "-":
            seq1_best += BLOSUM62[aa1, aa1]
        if aa2 != "-":
            seq2_best += BLOSUM62[aa2, aa2]

        if aa1 == "-" or aa2 == "-":
            if prev_aa1 == "-" or prev_aa2 == "-":
                observed_score += gap_extend
            else:
                observed_score += gap_open
        else:
            observed_score += BLOSUM62[aa1, aa2]
        prev_aa1 = str(aa1)
        prev_aa2 = str(aa2)

    subs_mat_score = ((observed_score / seq1_best) + (observed_score / seq1_best)) / 2

    # PSI PRED comparison
    num_gaps = 0
    ss_score = 0
    for row1 in psi_pred_files[id1].itertuples():
        if (psi_pred_files[id2]["indx"] == row1.indx).any():
            row2 = psi_pred_files[id2].loc[psi_pred_files[id2]["indx"] == row1.indx]
            row_score = 0
            row_score += 1 - abs(float(row1.coil_prob) - float(row2.coil_prob))
            row_score += 1 - abs(float(row1.helix_prob) - float(row2.helix_prob))
            row_score += 1 - abs(float(row1.sheet_prob) - float(row2.sheet_prob))
            ss_score += row_score / 3
        else:
            num_gaps += 1

    align_len = len(psi_pred_files[id2]) + num_gaps
    ss_score /= align_len
    final_score = (ss_score * 0.3) + (subs_mat_score * 0.7)
    with lock:
        with open(outfile, "a") as _ofile:
            _ofile.write("%s\t%s\t%s\n" % (id1, id2, final_score))
    return


def create_all_by_all_scores(seqs, group):
    printer.write("Running MAFFT")
    alignment = Alb.generate_msa(Sb.make_copy(seqs), tool="mafft", params="--globalpair --thread -1", quiet=True)
    printer.write("Updating PsiPred files")
    # Need to specify what columns the PsiPred files map to now that there are gaps.
    psi_pred_files = {}
    for rec in alignment.records_iter():
        ss_file = pd.read_csv("%s/psi_pred/%s.ss2" % (in_args.outdir, rec.id), comment="#",
                              header=None, delim_whitespace=True)
        ss_file.columns = ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]
        ss_counter = 0
        for indx, residue in enumerate(rec.seq):
            if residue != "-":
                ss_file.set_value(ss_counter, "indx", indx)
                ss_counter += 1
        psi_pred_files[rec.id] = ss_file

    printer.write("Removing gaps")
    alignment = Alb.trimal(alignment, "gappyout")
    # Re-update PsiPred files, now that some columns are removed
    for rec in alignment.records_iter():
        new_psi_pred = []
        for row in psi_pred_files[rec.id].itertuples():
            if alignment.alignments[0].position_map[int(row[1])][1]:
                new_psi_pred.append(list(row)[1:])
        psi_pred_files[rec.id] = pd.DataFrame(new_psi_pred, columns=["indx", "aa", "ss", "coil_prob",
                                                                     "helix_prob", "sheet_prob"])

    printer.write("Preparing to calculate Sim scores")
    alignment.write("%s/alignments/group_%s.aln" % (in_args.outdir, group))
    ids1 = [rec.id for rec in alignment.records_iter()]
    ids2 = [rec.id for rec in alignment.records_iter()]
    all_by_all = []
    for rec1 in ids1:
        del ids2[ids2.index(rec1)]
        for rec2 in ids2:
            all_by_all.append((rec1, rec2))
    outfile = "%s/sim_scores/group_%s.csv" % (in_args.outdir, group)
    printer.clear()
    MyFuncs.run_multicore_function(all_by_all, score_sequences, [alignment, psi_pred_files, outfile])
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
    parser.add_argument("-sr", "--supress_recursion", action="store_true",
                        help="Stop after a single round of MCL. For testing.")
    parser.add_argument("-scc", "--supress_clique_check", action="store_true",
                        help="Do not check for or break up cliques. For testing.")
    parser.add_argument("-ssf", "--supress_singlet_folding", action="store_true",
                        help="Do not check for or merge singlets. For testing.")
    parser.add_argument("-op", "--open_penalty", help="Penalty for opening a gap in pairwise alignment scoring",
                        type=float, default=-5)
    parser.add_argument("-ep", "--extend_penalty", help="Penalty for extending a gap in pairwise alignment scoring",
                        type=float, default=0)
    parser.add_argument("-nt", "--no_msa_trim", action="store_true",
                        help="Don't apply the gappyout algorithm to MSAs before scoring")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Overwrite previous run")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress all output during run (only final output is returned)")
    parser.add_argument("-psi", "--psi_pred", action="store",
                        help="Specify directory with psi_pred results")

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

    ambiguous_X = {"A": 0, "R": -1, "N": -1, "D": -1, "C": -2, "Q": -1, "E": -1, "G": -1, "H": -1, "I": -1, "L": -1,
                   "K": -1, "M": -1, "F": -1, "P": -2, "S": 0, "T": 0, "W": -2, "Y": -1, "V": -1}
    for aa in ambiguous_X:
        pair = sorted((aa, "X"))
        pair = tuple(pair)
        PHAT[pair] = ambiguous_X[aa]
        BLOSUM62[pair] = ambiguous_X[aa]
        BLOSUM45[pair] = ambiguous_X[aa]

    gap_open = in_args.open_penalty
    gap_extend = in_args.extend_penalty

    os.makedirs(in_args.outdir)
    os.makedirs("%s/alignments" % in_args.outdir)
    os.makedirs("%s/mcmcmc" % in_args.outdir)
    os.makedirs("%s/sim_scores" % in_args.outdir)
    os.makedirs("%s/psi_pred" % in_args.outdir)
    if in_args.psi_pred and os.path.isdir(in_args.psi_pred):
        files = os.listdir(in_args.psi_pred)
        for f in files:
            shutil.copyfile("%s/%s" % (in_args.psi_pred, f), "%spsi_pred/%s" % (in_args.outdir, f))

    print("\nExecuting PSI-Pred")
    MyFuncs.run_multicore_function(sequences.records, _psi_pred)

    print("\nGenerating initial all-by-all")
    scores_data = create_all_by_all_scores(sequences, group="0")

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
