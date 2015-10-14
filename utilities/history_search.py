#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Feb 24 2015 

"""
Return regex of my .history folder
"""

import os
import sys
import re
import argparse
from subprocess import Popen

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="history", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("regex", help="Pattern to search for in history", action="store", nargs='?', default=".*")
    parser.add_argument("depth", help="Number of matches to return", action="store", nargs='?', type=int, default=10)
    parser.add_argument("-d", "--date", help="Specify a specific history file", action="store")
    parser.add_argument("-r", "--run", help="Re-run the most recent result matching your search", action="store_true")
    parser.add_argument("-ld", "--list_dates", help="List the available history files", action="store_true")
    parser.add_argument("-pid", "--process_id", help="Restrict results to specific pid.", action="store")

    in_args = parser.parse_args()

    root, dirs, hist_files = next(os.walk("/Volumes/Zippy/.history/"))

    if in_args.list_dates:
        for _file in hist_files:
            print(_file)
        sys.exit()

    if in_args.date:
        hist_files = [_file for _file in hist_files if re.search(in_args.date, _file)]

    if not hist_files:
        sys.exit("No history files specified.")

    hist_files.sort()

    history_list = []
    for _file in hist_files:
        with open("%s/%s" % (root, _file), "r") as ifile:
            history_list += ifile.readlines()

    output = []
    while in_args.depth > 0 and len(history_list) > 0:
        line = history_list.pop()
        if in_args.process_id:
            if not re.search(":%s " % in_args.process_id, line):
                continue

        command = re.sub(":[0-9]* [JFMASOND][a-z]{2}/[0-3][0-9]/[0-9]{2} [0-9]{2}:[0-9]{2}; ", "", line)

        if re.search(in_args.regex, command):
            if in_args.run:
                print(">>> %s" % command)
                Popen(command, shell=True).wait()
                sys.exit()

            output.append(line)
            in_args.depth -= 1

    output.reverse()
    print("".join(output).strip() if output else "No history found for that search")
