#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Feb 4 2015 

"""
Simple program to do regular expression substitutions in a file.

To keep part of a match, surround the part you want to keep in parentheses, and then in the
replace string, include escaped integers for each part you want kept.

eg. ./regrep.py "<div name='([A-Za-z1-9]*)'>(.*)</div>" "<span id='\1'>\2</span>" infile.html
"""

import re
import argparse

parser = argparse.ArgumentParser(prog="regrep", description="Do a regular expression replace", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("regex", help="pattern to search", action="store")
parser.add_argument("replace", help="replacement", action="store")
parser.add_argument("file", help="input file", action="store")
parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')

in_args = parser.parse_args()

if __name__ == '__main__':
    with open(in_args.file, "r") as ifile:
        file_conts = ifile.read()

    sub_text = re.sub(in_args.regex, in_args.replace, file_conts, flags=re.MULTILINE)

    if in_args.in_place:
        with open(in_args.file, "w") as ofile:
            ofile.write(sub_text)

    else:
        print(sub_text)