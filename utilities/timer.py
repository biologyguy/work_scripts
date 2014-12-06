#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Nov 20 2014

"""
Set a timer
"""

from time import time, sleep
from MyFuncs import DynamicPrint, pretty_time
import argparse

parser = argparse.ArgumentParser(prog="timer.py", description="Simple command line timer",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('length', help='How many seconds to count?', type=int)
parser.add_argument('-q', '--quiet', help="Don't display the counter", action="store_true")

in_args = parser.parse_args()

start = round(time())
end = start + in_args.length
if not in_args.quiet:
    dynamic_print = DynamicPrint()
    while time() < end:
        dynamic_print.write(pretty_time(round(time()) - start))
        sleep(1)
    dynamic_print.write("%s" % pretty_time(in_args.length))

else:
    sleep(in_args.length)
