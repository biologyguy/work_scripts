#!/home/bondsr/bin/python3
# -*- coding: utf-8 -*-
# Created on: Oct 24 2014 

"""
Generates a bunch of trees for all partitions in an alignement. Wraps raxml.
"""

import argparse
import sys
import os
import shutil
import MyFuncs
from subprocess import Popen, PIPE


def ml_trees(aln_file):
    """Run standard RAxML on partitions"""
    part_name = aln_file.split(".")[2]
    command = "/home/bondsr/bin/raxmlHPC-MPI-SSE3 -s %s/%s -n %s_%s -m GTRGAMMAI -f d -p 12345 -N 20 -w %s/ML_trees/ --no-bfgs;" % (alignments, aln_file, base_name, part_name, work_dir)

    stdout = Popen(command, shell=True, stdout=PIPE).communicate()[0]

    with open("%s/ML_%s.log" % (qsub_files, part_name), "w") as ofile:
        ofile.write("Computing ML trees with RAxML\n%s\n" % command)
        ofile.write(stdout.decode("utf-8"))

    return


def bootstrap(aln_file):
    """Run standard RAxML on partitions"""
    part_name = aln_file.split(".")[2]
    command = "/home/bondsr/bin/raxmlHPC-MPI-SSE3 -s %s/%s -n %s_%s -m GTRGAMMAI -f d -p 12345 -b 12345 -N 100 -w %s/bootstraps/;" % (alignments, aln_file, base_name, part_name, work_dir)

    stdout = Popen(command, shell=True, stdout=PIPE).communicate()[0]

    with open("%s/ML_%s.log" % (qsub_files, part_name), "a") as ofile:
        ofile.write("\n\nComputing bootstraps with RAxML\n%s\n" % command)
        ofile.write(stdout.decode("utf-8"))

    return


def best_tree(aln_file):
    """Run standard RAxML on partitions"""
    part_name = aln_file.split(".")[2]
    command = "/home/bondsr/bin/best_ML_tree.py %s/ML_trees/RAxML_info.%s_%s;" % (work_dir, base_name, part_name)

    stdout = Popen(command, shell=True, stdout=PIPE).communicate()[0]

    with open("%s/ML_%s.log" % (qsub_files, part_name), "a") as ofile:
        ofile.write("\n\nFinding the best ML tree\n%s\n" % command)
        ofile.write(stdout.decode("utf-8"))

    return


def cons_tree(aln_file):
    """Run standard RAxML on partitions"""
    part_name = aln_file.split(".")[2]
    command = "/home/bondsr/bin/raxmlHPC-PTHREADS-SSE3 -n %s_%s -m GTRGAMMAI -f b -p 12345 -T 2 -w %s/cons_trees/ -t %s/ML_trees/RAxML_best.%s_%s -z %s/bootstraps/RAxML_bootstrap.%s_%s;" % (base_name, part_name, work_dir, work_dir, base_name, part_name, work_dir, base_name, part_name)

    stdout = Popen(command, shell=True, stdout=PIPE).communicate()[0]

    with open("%s/ML_%s.log" % (qsub_files, part_name), "a") as ofile:
        ofile.write("\n\nGenerating consensus trees\n%s\n" % command)
        ofile.write(stdout.decode("utf-8"))

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="raxml_qsub_pipeline", description="This script is currently very specific to my Gryphon home directory")

    parser.add_argument("partition_trees_dir", help="directory all the magic happens", action="store")
    #parser.add_argument("-t", "--true", help="", action="store_true", default=False)
    #parser.add_argument("-c", "--choice", help="", type=str, choices=["", ""], default=False)
    #parser.add_argument("-m", "--multi_arg", nargs="+", help="", default=[])

    in_args = parser.parse_args()

    work_dir = os.path.abspath(in_args.partition_trees_dir)
    os.chdir(work_dir)
    alignments = os.path.abspath("alignments/")
    ML_trees = os.path.abspath("ML_trees/")
    bootstraps = os.path.abspath("bootstraps/")
    cons_trees = os.path.abspath("cons_trees/")
    qsub_files = os.path.abspath("qsub_files/")
    outloop = False
    file_name = ""
    base_name = ""

    for (dirpath, dirnames, filenames) in MyFuncs.walklevel(os.path.abspath("../"), level=1):
        for file_name in filenames:
            if file_name.split(".")[-1] == "phy":
                outloop = True
                base_name = "_".join(file_name.split(".")[0].split("_")[1:])
                shutil.copy("%s/%s" % (os.path.abspath("../"), file_name), "alignments/%s.phy" % base_name)
                break

        if not outloop:
            sys.exit("ERROR: Failed to find the alignment file.")

        break

    command = "/home/bondsr/bin/raxmlHPC-PTHREADS-SSE3 -s %s/%s.phy -n %s -m GTRCAT -f s -T 2 -w %s/alignments/ -q %s/domain_partition.txt;" % (alignments, base_name, base_name, work_dir, work_dir)

    Popen(command, shell=True).wait()
    os.remove("%s/%s.phy" % (os.path.abspath("alignments/"), base_name))

    partition_list = []
    for (dirpath, dirnames, filenames) in MyFuncs.walklevel(os.path.abspath("alignments/"), level=1):
        for file_name in filenames:
            partition_list.append(file_name)
        break

    MyFuncs.run_multicore_function(partition_list, ml_trees)
    MyFuncs.run_multicore_function(partition_list, bootstrap)

    MyFuncs.run_multicore_function(partition_list, best_tree)
    MyFuncs.run_multicore_function(partition_list, cons_tree)
