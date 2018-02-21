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
import os
import sys
import argparse
from io import StringIO


def main():
    parser = argparse.ArgumentParser(prog="regrep", description="Do a regular expression replace", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input", help="input file", action="store", nargs="*", default=[sys.stdin])
    parser.add_argument("regex", help="pattern to search", action="store")
    parser.add_argument("replace", help="replacement", action="store")
    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')

    in_args = parser.parse_args()
    file_path = None
    in_args.input = in_args.input[0]

    if str(type(in_args.input)) == "<class '_io.TextIOWrapper'>":
        if not in_args.input.seekable():  # Deal with input streams (e.g., stdout pipes)
            input_txt = in_args.input.read()
            temp = StringIO(input_txt)
            in_args.input = temp
        in_args.input.seek(0)
        file_conts = in_args.input.read()

    elif type(in_args.input) == str and os.path.isfile(in_args.input):
        file_path = str(in_args.input)
        with open(file_path, "r") as ifile:
            file_conts = ifile.read()
    else:
        file_conts = in_args.input

    sub_text = re.sub(in_args.regex, in_args.replace, file_conts, flags=re.MULTILINE)

    if file_path and in_args.in_place:
        with open(file_path, "w") as ofile:
            ofile.write(sub_text)

    else:
        print(sub_text)


if __name__ == '__main__':
    main()
