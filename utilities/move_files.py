#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Oct 8 2017

"""
Search a folder and all sub-folders for files matching a regex, and then use flags to move or copy those files
somewhere else
"""

import sys
import os
import re
import shutil


# Pulled this function off of Stack Overflow -- posted by nosklo
# Iterates over directories only to a specified depth (useful in a for loop)
# Note that this is a generator, so need to use next() or `with` to get a result
def walklevel(some_dir, level):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if level != 0 and num_sep + level - 1 <= num_sep_this:
            del dirs[:]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="move_files", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("search_dir", action="store", help="Where do you want to search?")
    parser.add_argument("regex", action="append", nargs="+", help="One or more regular expressions to search for.")
    parser.add_argument("-c", "--copy", action="store", help="Specify a folder to copy matches to")
    parser.add_argument("-m", "--move", action="store", help="Specify a folder to move matches to")
    parser.add_argument("-l", "--level", type=int, action="store", default=0,
                        help="Restrict how deep to walk into sub-folders")
    in_args = parser.parse_args()

    if in_args.copy and in_args.move:
        sys.stderr.write("Error: Cannot 'copy' and 'move' at the same time. Pick one or the other")
        sys.exit()

    if (in_args.copy and not os.path.isdir(in_args.copy)) or (in_args.move and not os.path.isdir(in_args.move)):
        sys.stderr.write("Error: output directory does not exist. Please create it first.")
        sys.exit()

    files_to_move = []
    search_dir = os.path.abspath(in_args.search_dir)
    regex = "|".join(in_args.regex[0])
    for r, d, f in walklevel(search_dir, level=in_args.level):
        for _file in f:
            if re.search(regex, _file):
                files_to_move.append(os.path.join(r, _file))

    print("""
###################################
# Files identified in your search #
###################################
""")
    for _file in files_to_move:
        print(_file)

    if not in_args.copy and not in_args.move:
        print("\nUse the --move or --copy flag to actually modify your system")

    if in_args.copy:
        in_args.copy = os.path.abspath(in_args.copy)
        for _file in files_to_move:
            shutil.copy(_file, os.path.join(in_args.copy, os.path.split(_file)[1]))
        print("\nAll files copied to %s" % in_args.copy)

    if in_args.move:
        in_args.move = os.path.abspath(in_args.move)
        for _file in files_to_move:
            shutil.move(_file, os.path.join(in_args.move, os.path.split(_file)[1]))
        print("\nAll files moved to %s" % in_args.move)
