#!/usr/bin/env python3
# -*- coding: utf-8 -*-#
# Created on: Dec 7 2017

"""
Add a header to selected files
"""
import os
import re


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="add_header", description="Add a header to selected files",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("files", action="store", help="Regular expression to match files in current directory")
    parser.add_argument("header", action="store", help="String or path to file with text to be added to head of file")
    parser.add_argument("-e", "--escape_char", action="store", default="#", metavar="",
                        help="Specify the escape character to use")
    parser.add_argument("-c", "--commit", action="store_true", help="Actually modify the files")
    parser.add_argument("-s", "--simple", action="store_true", help="Only use a single instance of escape character")
    in_args = parser.parse_args()

    if os.path.isfile(os.path.abspath(in_args.header)):
        with open(os.path.abspath(in_args.header), "r") as ifile:
            in_args.header = [line.strip() for line in ifile.readlines()]
            print(in_args.header)
    else:
        in_args.header = in_args.header.split("\\n")

    max_line = max([len(line) for line in in_args.header])
    header = in_args.escape_char * (max_line + 4) + "\n" if not in_args.simple else ""
    for line in in_args.header:
        if in_args.simple:
            header += "{0} {1}\n".format(in_args.escape_char, line)
        else:
            header += "{0} {1} {0}\n".format(in_args.escape_char, line.ljust(max_line, " "))
    header += in_args.escape_char * (max_line + 4) if not in_args.simple else ""
    header = header.strip() + "\n"
    output = "\nHeader:\n\n" + header + "\n"

    mod_files = []
    root, dirs, files = next(os.walk(os.getcwd()))
    for next_file in files:
        if re.search(in_args.files, next_file):
            mod_files.append(next_file)

    if not mod_files:
        output += "No files were identified with the regular expression '%s'" % in_args.files
    elif in_args.commit:
        output += "Added to files:\n\n"
        for next_file in mod_files:
            with open(next_file, "r", encoding="utf8") as ifile:
                contents = ifile.read()
            with open(next_file, "w", encoding="utf8") as ofile:
                ofile.write("%s%s" % (header, contents))
    else:
        output += "Pass the `-c` flag to modify the following files:\n\n"

    if mod_files:
        output += "\t" + "\n\t".join(mod_files)
    output += "\n"
    print(output)
