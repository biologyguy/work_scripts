#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: May 30 2015 

"""
Convert the output from homolog_caller into an SVG tree of polytomies
"""
import sys

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
        output += "\ttree tree_1 = %s" % self.newick(self.node_list, 1)
        output += "\nend;"
        return output

    def newick(self, _nodes, depth):
        output = "("
        node_groups = {}
        for _node in _nodes:
            try:
                node_groups.setdefault(_node.rank[depth], []).append(_node)
            except IndexError:  # This happens when cliques are involved
                node_groups.setdefault(_node.rank[depth - 1], []).append(_node)

        node_groups = [_group for indx, _group in node_groups.items()]
        node_groups.sort(key=lambda x: x[0].support, reverse=True)

        for _group in node_groups:
            if len(_group) == 1:
                output += "("
                for _leaf in _group[0].leaf_list:
                    if not _leaf.support:
                        _leaf.support = 1.
                    output += "'%s':1.0," % _leaf.label
                output = "%s):1.0," % output.strip(",")

            else:
                output += self.newick(_group, depth + 1)

        output = "%s):1.0," % output.strip(",")
        if depth == 1:
            output = "%s;" % output.strip(",")
        return output

    def __len__(self):
        length = len([_leaf.label for _node in self.node_list for _leaf in _node.leaf_list])
        return length


class Node:
    def __init__(self, leaf_list, rank, ave_support, std_support):
        self.leaf_list = sorted(leaf_list, key=lambda x: x.support, reverse=True)
        self.rank = rank.split("_")
        self.support = ave_support
        self.std = std_support

    def __len__(self):
        return len(self.leaf_list)


class Leaf:
    def __init__(self, label, _support=None):
        self.label = label
        if _support == "nan":
            _support = None
        self.support = _support
        self.color = support_color(_support)


def support_color(_support):
    if not _support:
        color = "000000"

    else:
        _support = float(_support)
        if 1 < _support < 0:
            raise ValueError("Leaf support values must be between 0 and 1")

        red = hex_string(175 - (175 * _support)) if _support <= 0.5 else hex_string(0)
        green = hex_string(175 * _support) if _support >= 0.5 else hex_string(0)
        blue = hex_string(175 * (1 - (abs(0.5 - _support))))
        color = "%s%s%s" % (red, green, blue)

    return color


def hex_string(value):
    output = "0x%0.2X" % int(value)
    return output[2:]

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="homolog_tree_builder", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("support", help="", action="store")

    in_args = parser.parse_args()

    with open(in_args.support, "r") as ifile:
        support_file = ifile.read().strip().split("group_")[1:]

    nodes = []

    for node in support_file:
        node = node.strip().split("\t")
        support = node[1]
        leaves = []
        for leaf in node[2:]:
            leaves.append(Leaf(leaf, 1.0))

        nodes.append(Node(leaves, rank=node[0], ave_support=0,
                          std_support=0))

    nexus = Nexus(nodes)
    print(nexus.print())
