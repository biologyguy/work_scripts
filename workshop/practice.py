#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Sep 19 2016

# input file
InFileName = "/Users/jonesalm/practice_data/tps1.fa"

# open file in read mode
InFile = open(InFileName, 'r')
# create empty dictionary
SeqDict = {}

for Line in InFile:
    if Line.startswith('>'):
        # if the line starts with a '>', strip the line endings
        Line = Line.strip('\n')
        # split the line whenever there is a '|'
        ElementList = Line.split('|')
        # 4th element in the list is the access number
        AccessNum = ElementList[3]
        if AccessNum.startswith("NP"):
            # create a key in the SeqDict for each access number found
            SeqDict[AccessNum] = []
            print(">" + AccessNum)
    else:
        if AccessNum.startswith("NP"):
            # if the line doesnt start with '>', AKA a line of protein sequences, strip the line ending
            SeqLine = Line.strip('\n')
            # append this SeqLine to a list of seqs for that AccessNum
            SeqDict[AccessNum].append(SeqLine)
            print(SeqLine)

DictLength = len(SeqDict)
print(DictLength)
print('Number of Sequences: %d' % DictLength)

InFile.close()
