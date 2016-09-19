#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Sep 19 2016


def protein_parse(infile):

    # open file in read mode
    infile = open(infile, 'r')
    # create empty dictionary
    seq_dict = {}

    for line in infile:
        if line.startswith('>'):
            # if the line starts with a '>', strip the line endings
            line = line.strip('\n')
            # split the line whenever there is a '|'
            element_list = line.split('|')
            # 4th element in the list is the access number
            access_num = element_list[3]
            if access_num.startswith("NP"):
                # create a key in the seq_dict for each access number found
                seq_dict[access_num] = []
                print(">" + access_num)
        else:
            if access_num.startswith("NP"):
                # if the line doesnt start with '>', AKA a line of protein sequences, strip the line ending
                seq_line = line.strip('\n')
                # append this seq_line to a list of seqs for that access_num
                seq_dict[access_num].append(seq_line)
                print(seq_line)

    dict_length = len(seq_dict)
    print(dict_length)
    print('Number of Sequences: %d' % dict_length)

    infile.close()

if __name__ == '__main__':
    protein_parse('/Users/jonesalm/practice_data/tps1.fa')