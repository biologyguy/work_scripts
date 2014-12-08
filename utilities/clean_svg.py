#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Oct 23 2014 

"""
SVG files made by a lot of programs include font-family definitions that Illustrator chokes on. Uses this program to
strip them out.
"""

import argparse
import os
import re
import html

parser = argparse.ArgumentParser(prog='clean_svg', description='Fix SVG files so they will open in Illustrator')

parser.add_argument('svg_file', help='Location of the file you want to clean', action='store')
parser.add_argument('-o', '--outfile', help="Save the new file somewhere, otherwise it's sent to stdout.", action='store', default=False)
parser.add_argument('-i', '--in_place', help="Over-write original file", action='store_true')

in_args = parser.parse_args()
in_file = os.path.abspath(in_args.svg_file)
with open(in_file, "r") as ifile:
    content = ifile.read()

content = html.unescape(content)
content = re.sub(" ?font-family:.*?;", "", content)
content = re.sub("Œ±", "α", content)
content = re.sub("Œ≤", "β", content)
content = re.sub("Œ≥", "γ", content)

if in_args.in_place:
    with open(in_file, "w") as ofile:
        ofile.write(content)

elif in_args.outfile:
    with open(os.path.abspath(in_args.outfile), "w") as ofile:
        ofile.write(content)

else:
    print(content)