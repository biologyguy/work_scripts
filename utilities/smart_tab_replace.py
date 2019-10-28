#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Oct 6 2019

"""
Simple program to swap tab characters for the correct number of spaces (assuming tabs = 4 spaces).

eg. smart_tab_replace.py infile.java
"""

import re
import os
import sys
import argparse
import MyFuncs as mf


def main():
    parser = argparse.ArgumentParser(prog="smart_tab_replace", description="Clean up tabs with up to 4 spaces",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input", help="input file", action="store", nargs="*", default=[sys.stdin])
    parser.add_argument("-i", "--in_place", help="Rewrite the input file in-place. Be careful!", action='store_true')

    in_args = parser.parse_args()

    text = mf.easy_read(in_args.input[0])
    text = text.split("\n")

    for lindx, line in enumerate(text):
        tabs = list(re.finditer(r'\t', line))
        new_line = line if not tabs else line[:tabs[0].start()]
        for tindx, t in enumerate(tabs[1:]):
            tindx += 1
            new_line += " " * (4 - (len(new_line) % 4))
            new_line += line[tabs[tindx - 1].end():t.start()]
        if tabs:
            new_line += " " * (4 - (len(new_line) % 4))
            new_line += line[tabs[-1].end():]
        text[lindx] = new_line

    text = "\n".join(text)

    # Strip trailing whitespace
    text = re.sub(r' +$', '', text, flags=re.MULTILINE)

    # Remove excessive line breaks
    text = re.sub(r'\n\n+', r'\n\n', text)

    # Strip EOF line breaks
    text = re.sub(r'\n*$', '\n', text)

    file_path = None
    if type(in_args.input[0]) == str and os.path.isfile(in_args.input[0]):
        file_path = str(in_args.input[0])

    if file_path and in_args.in_place:
        with open(file_path, "w") as ofile:
            ofile.write(text)
    else:
        sys.stdout.write(text)


if __name__ == '__main__':
    main()
