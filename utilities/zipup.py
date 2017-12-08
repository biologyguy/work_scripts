#!/usr/bin/env python3
# -*- coding: utf-8 -*-#
# Created on: Dec 7 2017

"""
Make a new zip file with the included files
"""
import os
import re
from zipfile import ZipFile


def get_all_files(root):
    root, dirs, files = next(os.walk(os.path.abspath(root)))
    paths = [os.path.join(root, f) for f in files]
    for d in dirs:
        paths += get_all_files(os.path.join(root, d))
    return paths


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="zipup", description="Make a new zip file with the included files",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("files", action="append", nargs="+",
                        help="Regular expression(s) to match files in current directory and sub-directories")
    parser.add_argument("-o", "--outfile", action="store", default="new_zipup.zip", metavar="",
                        help="Name of new zip file")
    parser.add_argument("-c", "--commit", action="store_true", help="Actually modify the files")
    in_args = parser.parse_args()

    included = []
    cwd_size = len(os.getcwd() + os.sep)
    for next_file in get_all_files(os.getcwd()):
        if re.search("|".join(in_args.files[0]), os.path.split(next_file)[1]):
            included.append(next_file[cwd_size:])

    if in_args.commit:
        try:
            with ZipFile(in_args.outfile, mode="x") as zipfile:
                [zipfile.write(f) for f in included]
            output = "The following file(s) have been added to %s:\n\n\t%s\n" % (in_args.outfile, "\n\t".join(included))

        except FileExistsError:
            output = "Error: The file you are trying to create already exists. Pick a different name."

    else:
        output = "The following file(s) will be added to %s if you " \
                 "include the -c flag:\n\n\t%s\n" % (in_args.outfile, "\n\t".join(included))

    print(output)
