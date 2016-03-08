#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Feb 29 2016 

"""
Load up a bucket of data into RAM, then sit on it until you're finished stressing your computer.
"""
from time import sleep
import psutil

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="stress", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Load up a bucket of data into RAM, then sit on it until you're "
                                                 "finished stressing your computer.", )

    parser.add_argument("-m", "--max", type=int, action="store", help="Max GB to fill up.")

    in_args = parser.parse_args()

    if not in_args.max:
        memory = psutil.virtual_memory()
        in_args.max = int(round(memory[0] / 1000000))

    else:
        in_args.max *= 1000

    megabyte = (0,) * int((1024 * 1024 / 8))
    print("Loading %s GB of data into RAM" % (in_args.max / 1000))
    data = megabyte * in_args.max

    print("Going to sleep now...")
    try:
        while True:
            sleep(10)

    except KeyboardInterrupt:
        print("Freeing up memory.")
        data = 1

    print("Goodbye")
