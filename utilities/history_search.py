#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Feb 24 2015

"""
Return regex of my .history folder
"""

import argparse
import os
import re
from subprocess import Popen, PIPE
import sys
import time

if __name__ == '__main__':
    rows, columns = os.popen('stty size', 'r').read().split()  # Get terminal window size

    parser = argparse.ArgumentParser(prog="history", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("regex", help="Pattern to search for in history", action="store", nargs='?', default=".*")
    parser.add_argument("depth", help="Number of matches to return", action="store", nargs='?', type=int, default=10)
    parser.add_argument("-a", "--all", help="Show all results from given search.", action="store_true")
    parser.add_argument("-d", "--date", help="Specify a specific history file", action="store")
    parser.add_argument("-dlt", "--delete", help="Remove history records", type=int, action="store")
    parser.add_argument("-e", "--expand", help="Fully expand all returned commands", action="store_true")
    parser.add_argument("-hist", "--history", help="Look at hist calls", action="store_true")
    parser.add_argument("-l", "--length", default=columns, type=int,
                        help="Length of commands returned. Default is set to width of terminal.")
    parser.add_argument("-r", "--run", nargs="?", type=int, action="append",
                        help="Re-run the most recent result matching your search")
    parser.add_argument("-rc", "--run_careful", nargs="?", type=int, action="append",
                        help="Re-run the most recent result matching your search, but check first")
    parser.add_argument("-u", "--unique", action="store_true", help="Do not show replicate commands")
    parser.add_argument("-ld", "--list_dates", help="List the available history files", action="store_true")
    parser.add_argument("-pid", "--process_id", help="Restrict results to specific pid.", action="store")

    in_args = parser.parse_args()

    save_dir = Popen("echo $HISTORYREC", shell=True, stdout=PIPE).communicate()[0].decode().strip()
    root, dirs, hist_files = next(os.walk("%s/.history/" % save_dir))

    careful = False
    if in_args.run_careful:
        careful = True
        in_args.run = in_args.run_careful

    if in_args.run:
        in_args.run = 1 if not in_args.run[0] else in_args.run[0]

    if in_args.history:
        hist_files = [_file for _file in hist_files if _file.startswith("hist")]
    else:
        hist_files = [_file for _file in hist_files if _file.startswith("bash")]

    if in_args.list_dates:
        for _file in hist_files:
            print(_file)
        sys.exit()

    if in_args.date:
        hist_files = [_file for _file in hist_files if re.search(in_args.date, _file)]

    if not hist_files:
        sys.exit("No history files specified.")

    hist_files.sort(reverse=True)

    output = []
    commands = []
    counter = 1
    breakout = False

    for _file in hist_files:
        line_num = -1
        if breakout:
            break

        with open("%s/%s" % (root, _file), "r") as ifile:
            history_list = ifile.readlines()

        while len(history_list) > 0:
            line = history_list.pop()
            line_num += 1
            if in_args.process_id:
                if not re.search(":%s " % in_args.process_id, line):
                    continue

            command = re.sub(":[0-9]* [JFMASOND][a-z]{2}/[0-3][0-9]/[0-9]{2} [0-9]{2}:[0-9]{2}; ", "", line)
            if in_args.unique and command in commands:
                continue
            commands.append(command)

            if re.search(in_args.regex, command):
                if in_args.run and in_args.run == counter:
                    print(">>> %s" % command)
                    if not careful or input("Continue [no]? ").lower().startswith("y"):
                        try:
                            Popen(command, shell=True).wait()
                        except KeyboardInterrupt:
                            pass

                    # Write separate history files for cases where history_search is used to run a command
                    date = time.strftime("%b/%d/%y %H:%M")
                    output = ":%s %s; %s" % (os.getppid(), date, command)
                    current_file = "hist_history-%s" % time.strftime("%Y%m")

                    with open("%s/.history/%s" % (save_dir, current_file), "a") as ofile:
                        ofile.write("%s\n" % output.strip())

                    sys.exit()

                if in_args.delete and in_args.delete == counter:
                    with open("%s/%s" % (root, _file), "r") as ifile:
                        lines = ifile.readlines()
                        lines.reverse()
                        lines = lines[:line_num] + lines[line_num+1:]
                        lines.reverse()
                    with open("%s/%s" % (root, _file), "w") as ofile:
                        print(">>> DELETING: %s" % line)
                        ofile.write("".join(lines))
                    sys.exit()

                line = str(counter) + line

                if len(line) > in_args.length and not in_args.expand:
                    left = round((in_args.length - 5) * 0.75)
                    right = round((in_args.length - 5) * 0.25) * -1
                    line = "%s ... %s" % (line[:left], line[right:])

                output.append(line)
                if not in_args.all:
                    in_args.depth -= 1
                counter += 1

                if in_args.depth < 1:
                    breakout = True
                    break

    output.reverse()
    print("".join(output).strip() if output else "No history found for that search")
