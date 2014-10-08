#!/usr/bin/python3
# -*- coding: utf-8 -*-

from Bio import AlignIO
import argparse, sys

parser = argparse.ArgumentParser(prog="Make partitionfinder cfg file", description="",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', help='Location of input phylip file', action='store')
parser.add_argument("-b", "--blocks", nargs='+', help='space separated list of block positions (e.g., 1-982 983-2193 3453-3876). If not specified, the entire sequence is used.', default=False)
in_args = parser.parse_args()

in_file = in_args.input.split("/")[-1]
gene_ids = []

with open(in_args.input, "r") as ifile:
    alignment = next(AlignIO.parse(ifile, "phylip-relaxed"))
    sequences = alignment.get_all_seqs()

    for seq in sequences:
        gene_ids.append(seq.id)

output = "## ALIGNMENT FILE ##\n"
output += "alignment = %s;\n\n" % in_file

output += "## BRANCHLENGTHS: linked | unlinked ##\n"
output += "branchlengths = linked;\n\n"  # options linked/unlinked

output += "## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | beast | <list> ##\n"
output += "##              for PartitionFinderProtein: all_protein | <list> ##\n"
output += "models = all;\n\n"  # options all/some group pf models to test

output += "# MODEL SELECCTION: AIC | AICc | BIC #\n"
output += "model_selection = AIC;\n\n"

output += "## DATA BLOCKS: see manual for how to define ##\n"
output += "[data_blocks]\n"

if in_args.blocks:
    blocks = in_args.blocks
else:
    blocks = ["1-%s" % alignment.get_alignment_length()]

for block in blocks:
    for gene in gene_ids:
        output += "%s = %s\\3;\n" % (gene, block)


output += "\n## SCHEMES, search: all | greedy | rcluster | hcluster | user ##\n"
output += "[schemes]\n"
output += "search = greedy;\n\n"

#output += "#user schemes go here if search=%s#" % ??

print(output)