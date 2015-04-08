#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import os
from Bio import SeqIO
from scipy import stats
import numpy as np
from math import ceil, log10
from run_multicore_function import run_multicore_function
from multiprocessing import Lock


class shuffled():

    def __init__(self, out_dir_path):
        self.prob_list = []

        if out_dir_path[-1:] == "/":
            out_dir_path = out_dir_path[0:-1]

        self.out_dir_path = out_dir_path

        trim_files = []
        for files in os.listdir(self.out_dir_path):
            if not re.search("_trim", files):
                trim_files.append(files)

        if len(trim_files) != 2:
            sys.exit(
                "Couldn't find exactly 2 _trim.fasta files in your output dir. Did you point me in the right direction?")

        fasta_file1 = "%s/%s" % (out_dir_path, trim_files[0])
        fasta_file2 = "%s/%s" % (out_dir_path, trim_files[1])

        with open(fasta_file1, "r") as fasta_file:
            self.fasta1_list = list(SeqIO.parse(fasta_file, 'fasta'))

        with open(fasta_file2, "r") as fasta_file:
            self.fasta2_list = list(SeqIO.parse(fasta_file, 'fasta'))

        # determine naming order
        for files in os.listdir(self.out_dir_path + "/ALIGNME_FILES/"):
            if not re.match(self.fasta2_list[0].id, files):
                tmp = self.fasta1_list
                self.fasta1_list = self.fasta2_list
                self.fasta2_list = tmp
                break

    def print_probabilities(self):
        # Having looked at the distributions of large sets of randomly shuffled alignment scores, the data is multimodal
        # Estimate the Gaussian-kernel-density, and then integrate from 0 to the un-shuffled score to get the probability
        # of observing a score at least as small as the one seen.
        print "Compiling shuffle scores and P-values"

        for rec1 in self.fasta1_list:
            for rec2 in self.fasta2_list:
                prob = self._get_pairwise_P_value(rec1.id, rec2.id)
                self.prob_list.append(prob["prob_0_to_score"])

        with open("%s/probabilities.dat" % self.out_dir_path, "w") as file:
            for _next in self.prob_list:
                neg_log = (-1) * log10(_next)
                file.write("%s\t%s\n" % (_next, neg_log))

        box_plot = []
        return True

    def data_output(self, data_list, split_out_file, comb_out_file):
        # output all bootstrap data for each seq-seq comparison
        # output a single column of combined data, using shuffled score minus
        # unshuffled score for all samples
        comb_out_file = open(comb_out_file, "w")
        combined_data = []

        split_out_file = open(split_out_file, "w")
        split_out_file.write("#")

        for next_pair in data_list:
            split_out_file.write(
                next_pair[0] + "_Orig\t" + next_pair[0] + "_Boots\t")
        split_out_file.write("\n")

        for next_pair in data_list:
            split_out_file.write(next_pair[1] + "\t" + next_pair[2][0] + "\t")
            comb_out_file.write(
                str(float(next_pair[1]) - float(next_pair[2][0])) + "\n")
            combined_data.append(
                str(float(next_pair[1]) - float(next_pair[2][0])))

        split_out_file.write("\n")

        for index in range(99):
            for next_pair in data_list:
                split_out_file.write("\t" + next_pair[2][index + 1] + "\t")
                comb_out_file.write(
                    str(float(next_pair[1]) - float(next_pair[2][index + 1])) + "\n")
                combined_data.append(
                    str(float(next_pair[1]) - float(next_pair[2][index + 1])))
            split_out_file.write("\n")
        split_out_file.close()
        comb_out_file.close()

        return combined_data

    def _percentile(self, val_list, perc):
        val_list = [float(i) for i in val_list]
        val_list = sorted(val_list)
        index = ceil(float(len(val_list)) * perc)
        return val_list[int(index)]

    def _score_alignme(self, alignme_file):
        file_lines = alignme_file.readlines()

        # clear out header rows
        while True:
            if re.match("#", file_lines[0]) != None:
                del file_lines[0]
            else:
                break
        del file_lines[-1]

        tally = 0.0
        count = 0
        for next in file_lines:
            regular = re.sub("\s+", ",", next)
            regular = re.sub("\?0", "0", regular)
            data = regular.split(",")
            tally += abs(float(data[1]) - float(data[5])) + abs(float(data[2]) - float(data[6])) + abs(
                float(data[3]) - float(data[7])) + abs(float(data[4]) - float(data[8]))

        return round(tally, 1)

    def _get_pairwise_P_value(self, id1, id2):
        output = {}
        output["pair"] = "%s-%s" % (id1, id2)

        with open("%s/ALIGNME_FILES/%s-%s.prf" % (self.out_dir_path, id1, id2), "r") as initial:
            output["score"] = self._score_alignme(initial)

        output["shuffled_scores"] = []
        output["shuff_minus_unshuff"] = []
        for next in range(100):
            file = "%s/ALIGNME_FILES/SHUFFLED_FILES/%s-%s/%s_%04d-%s_%04d.prf" % (
                self.out_dir_path, id1, id2, id1, next, id2, next)
            with open(file, "r") as prf:
                shuf_score = self._score_alignme(prf)
                output["shuffled_scores"].append(shuf_score)
                output["shuff_minus_unshuff"].append(
                    round(shuf_score - output["score"], 1))

        # check for normal
        output["shapiro"] = stats.shapiro(output["shuffled_scores"])

        # Build probability density function to describe the sampe data, and integrate the area under the curve for
        # ≤ unshuffled score
        kernal = stats.kde.gaussian_kde(output["shuffled_scores"])
        output["prob_0_to_score"] = kernal.integrate_box_1d(0, output["score"])
        return output

    '''
    #combined_data = data_output(output, self.out_dir_path + "/ALIGNME_FILES/bootstraps.dat",self.out_dir_path + "/ALIGNME_FILES/shuf_minus_unshuf.dat")
    #output data as self._percentiles
    out_file = open(self.out_dir_path + "/ALIGNME_FILES/bootstraps.perc","w")
    out_file.write("#Num\tGroup\t10th\t25th\t50th\t75th\t90th\tUnshuffled\n")
    counter = 1

    for next_pair in output:
        out_file.write(str(counter) + "\t" + next_pair[0] + "\t" + str(self._percentile(next_pair[2],0.1)) + "\t" + str(self._percentile(next_pair[2],0.25)) + "\t" + str(self._percentile(next_pair[2],0.5)) + "\t" + str(self._percentile(next_pair[2],0.75)) + "\t" + str(self._percentile(next_pair[2],0.9)) + "\t" + next_pair[1] + "\n")
        counter += 1
    out_file.close()

    out_file = open(self.out_dir_path + "/ALIGNME_FILES/combined_bootstraps.perc","w")
    out_file.write("#Num\tGroup\t10th\t25th\t50th\t75th\t90th\n")

    out_file.write(str(1) + "\tEverything\t" + str(self._percentile(combined_data,0.1)) + "\t" + str(self._percentile(combined_data,0.25)) + "\t" + str(self._percentile(combined_data,0.5)) + "\t" + str(self._percentile(combined_data,0.75)) + "\t" + str(self._percentile(combined_data,0.9)))
    out_file.close()
    '''