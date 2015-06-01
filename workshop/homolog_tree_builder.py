#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: May 30 2015 

"""
Convert the output from homolog_caller into an SVG tree of polytomies
"""

import sys
import os
import re
import shutil
import MyFuncs
import SeqBuddy
import svgwrite

class Nexus:
    def __init__(self, node_list):
        self.node_list = node_list

    def print(self):
        output = "#NEXUS\n"
        output += "begin taxa;\n"
        output += "\tdimensions ntax=%s;\n" % len(self)
        output += "\ttaxlabels\n"
        for _node in self.node_list:
            for _leaf in _node.leaf_list:
                output += "\t'%s'[&!color=#%s]\n" % (_leaf.label, _leaf.color)
        output += ";\n"
        output += "end;\n\n"
        output += "begin trees;\n"
        output += "\ttree tree_1 = %s" % self.newick()

        return output

    def newick(self):
        output = ""
        test = []
        for _node in self.node_list:
            output += "("
            test.append(_node.rank)
            for _leaf in _node.leaf_list:
                output += "%s %s,"
        print(test)
        test.sort()
        sys.exit(test)
        return output

    def __len__(self):
        length = len([leaf.label for node in self.node_list for leaf in node.leaf_list])
        return length

class Node:
    def __init__(self, leaf_list, rank, ave_support, std_support):
        self.leaf_list = leaf_list
        self.rank = rank
        self.support = ave_support
        self.std = std_support

    def __len__(self):
        return len(self.leaf_list)

class Leaf:
    def __init__(self, label, _support=None):
        self.label = label
        if not _support or _support == "nan":
            self.color = "000000"
            self.support = None

        else:
            self.support = float(_support)
            if 1 < self.support < 0:
                raise ValueError("Leaf support values must be between 0 and 1")

            red = hex_string(255 - (255 * self.support))
            green = hex_string(255 * self.support)
            self.color = "%s%s00" % (red, green)


def hex_string(value):
    output = "0x%0.2X" % value
    return output[2:]

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="homolog_tree_builder", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("support", help="", action="store")
    #parser.add_argument("-t", "--true", help="", action="store_true", default=False)
    #parser.add_argument("-c", "--choice", help="", type=str, choices=["", ""], default=False)
    #parser.add_argument("-m", "--multi_arg", nargs="+", help="", default=[])

    in_args = parser.parse_args()

    with open(in_args.support, "r") as ifile:
        support_file = ifile.read().strip().split("group_")[1:]

    nodes = []
    for node in support_file:
        node = node.strip().split("\n")
        support = node.pop(0).split(" ")
        leaves = []
        for leaf in node:
            leaf = leaf.split(" ")
            leaves.append(Leaf(leaf[0], leaf[1]))

        nodes.append(Node(leaves, rank=support[0], ave_support=float(support[1]),
                          std_support=float(support[2])))

    nexus = Nexus(nodes)
    print(nexus.print())

    """
    dwg = svgwrite.Drawing('test_files/test.svg', profile='tiny')
    dwg.add(dwg.line(start=(10, 10), end=(100, 100), stroke=svgwrite.rgb(10, 10, 16, '%')))
    dwg.add(dwg.line(start=(20, 20), end=(120, 120), stroke=svgwrite.rgb(10, 10, 16, '%')))
    dwg.add(dwg.text('Test', insert=(20, 20)))
    dwg.save()
    """