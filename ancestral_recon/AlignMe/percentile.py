#!/usr/bin/python
# -*- coding: utf-8 -*- 

import argparse, sys, random
import numpy as np
from math import ceil
from scipy import stats

def percentile(val_list,perc):
    val_list = [float(i) for i in val_list]
    val_list = sorted(val_list)
    index = ceil(float(len(val_list))*perc)
    return(val_list[int(index)])




parser = argparse.ArgumentParser(prog="", description="");
parser.add_argument('in_file', help='', action='store');
parser.add_argument('output_dir', help='', action='store');
#parser.add_argument('perc',help='percentile, from 0.0-1.0')
in_args = parser.parse_args()

if in_args.output_dir[-1:] == "/":
    in_args.output_dir = in_args.output_dir[0:-1]

with open(in_args.in_file,"r") as file:
    data = file.readlines()


out_file = open(in_args.output_dir + "/ALIGNME_FILES/combined_bootstraps.perc","w")
out_file.write("#Num\tGroup\t10th\t25th\t50th\t75th\t90th\n")

out_file.write(str(1) + "\tEverything\t" + str(percentile(data,0.1)) + "\t" + str(percentile(data,0.25)) + "\t" + str(percentile(data,0.5)) + "\t" + str(percentile(data,0.75)) + "\t" + str(percentile(data,0.9)))    
out_file.close()


data = [float(i) for i in data]
print "Average: " + str(round(np.mean(data,axis=0),4))
print "StDev: " + str(round(np.std(data,axis=0),4))
print "1-Samp T-test: %.3f. p-value: %.3g." % stats.ttest_1samp(data,0.0)
print "Normal test. Z-score: %.3g, P-value: %.3g " % stats.normaltest(data)

random.shuffle(data)
print "1-Samp T-test: %.3f. p-value: %.3g." % stats.ttest_1samp(data[0:1000],0.0)

random.shuffle(data)
print "1-Samp T-test: %.3f. p-value: %.3g." % stats.ttest_1samp(data[0:1000],0.0)

random.shuffle(data)
print "1-Samp T-test: %.3f. p-value: %.3g." % stats.ttest_1samp(data[0:1000],0.0)