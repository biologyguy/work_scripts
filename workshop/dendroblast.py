
import argparse
from scipy.cluster.hierarchy import linkage, to_tree, dendrogram
import numpy as np
import matplotlib.pylab as plt
import SeqBuddy
import sys

parser = argparse.ArgumentParser(prog="dendroblast.py", description="")
parser.add_argument("sequence", help="Supply a file path or a raw sequence", nargs="+", default=sys.stdin)
parser.add_argument("-drb", "--dendroblast", action="store_true",
                    help="Create a dendrogram from pairwise blast bit-scores. Returns newick format.")

in_args = parser.parse_args()


def denroblast(_seqs):  # This does not work yet... See Kelly and Maini, 2013, PlosONE
    _blast_res = SeqBuddy.bl2seq(_seqs).split("\n")

    # format the data into a dictionary for easier manipulation
    dist_dict = {}
    for pair in _blast_res:
        pair = pair.split("\t")
        if pair[0] in dist_dict:
            dist_dict[pair[0]][pair[1]] = float(pair[5])
        else:
            dist_dict[pair[0]] = {pair[1]: float(pair[5])}

        if pair[1] in dist_dict:
            dist_dict[pair[1]][pair[0]] = float(pair[5])
        else:
            dist_dict[pair[1]] = {pair[0]: float(pair[5])}

    # make a numpy array and headings list
    dist_array = np.zeros([len(_seqs.seqs), len(_seqs.seqs)])
    headings = []
    i = 0
    for _seq1 in _seqs.seqs:
        headings.append(_seq1.id)
        j = 0
        for _seq2 in _seqs.seqs:
            if _seq1.id != _seq2.id:
                dist_array[i][j] = dist_dict[_seq1.id][_seq2.id]
            j += 1
        i += 1

    headings.reverse()
    print(headings)
    data_link = linkage(dist_array, method='complete')
    dendrogram(data_link, labels=headings)
    plt.xticks(fontsize=8, rotation=90)
    plt.savefig("dendrogram.svg", format='svg')
    plt.show()

sequence = SeqBuddy.SeqBuddy(in_args.sequence)

denroblast(sequence)