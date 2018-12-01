#!/usr/bin/env python
from subprocess import Popen, PIPE
from os.path import join
from buddysuite import buddy_resources as br
import os
import string
from random import choice
import shutil


def convert(f, args):
    in_place = args[0]
    tmp = br.TempDir()
    rand_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
    tmp_file = join(tmp.path, rand_name + ".flac")
    out_file = join(tmp.path, rand_name + ".mp3")
    shutil.copyfile(f, tmp_file)

    Popen("ffmpeg -i '%s' '%s'" % (tmp_file, out_file), stderr=PIPE, stdout=PIPE, shell=True).communicate()
    shutil.copyfile(out_file, os.path.splitext(f)[0] + ".mp3")
    if in_place:
        os.remove(f)
    return


def main():
    import argparse
    parser = argparse.ArgumentParser(prog="flac2mp3", description="Wrapper for ffmpeg",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("music_dirs", action="append", nargs='+', help="Where are the files?")
    parser.add_argument("-i", "--in_place", action="store_true", help="Replace the flac files")

    in_args = parser.parse_args()

    music_dirs = [os.path.abspath(d) for d in in_args.music_dirs[0]]
    flac_files = []
    for d in music_dirs:
        root, dirs, files = next(os.walk(d))
        flac_files += [join(root, f) for f in files if os.path.splitext(f)[1].lower() == ".flac"]

    br.run_multicore_function(flac_files, convert, func_args=[in_args.in_place])


if __name__ == '__main__':
    main()
