#!/usr/bin/env python

from subprocess import Popen
import argparse
import os
import shutil


def encrypt(infile):
    Popen("gpg -e -u biologyguy@gmail.com -r biologyguy@gmail.com %s" % infile, shell=True).communicate()
    shutil.move(infile + ".gpg", infile)


def main():
    parser = argparse.ArgumentParser(prog="decrypt", description="Decrypt a file")

    parser.add_argument("infile", help="input file", action="store")
    parser.add_argument("-v", "--vim", help="Modify the file in vim, then encrypt again", action="store_true")
    parser.add_argument("-p", "--print", help="Print the contents of the file, then encrypt again", action='store_true')

    in_args = parser.parse_args()

    Popen("gpg -o {0}.gpg -d {0}".format(in_args.infile), shell=True).wait()
    if os.path.isfile(in_args.infile + ".gpg"):
        shutil.move(in_args.infile + ".gpg", in_args.infile)

    if in_args.vim:
        Popen("vim %s" % in_args.infile, shell=True).wait()
        encrypt(in_args.infile)

    elif in_args.print:
        with open(in_args.infile, "r") as ifile:
            print(ifile.read())
        encrypt(in_args.infile)


if __name__ == '__main__':
    main()
