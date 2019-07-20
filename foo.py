#!/usr/bin/env python3

import argparse
import utilities as br
import re
import sys
from math import floor
from random import randint


class Bug(object):
    def __init__(self, rgb):
        self.rgb = rgb
        self.r = rgb[0]
        self.g = rgb[1]
        self.b = rgb[2]

        self.energy = 100


def fight(bug1, bug2):
    bug1_power = 0
    # Red vs. Green (Red is always 10% stronger)
    RvG1 = (((bug1.r - ((bug1.g + bug1.b) / 10)) / 255) * bug2.g) * 1.1
    bug1_power += RvG1
    print(round(bug1_power, 2))

    # Green vs. Blue (Blue is always 10% more defensive)
    GvB1 = (((bug1.g - ((bug1.b + bug1.r) / 10)) / 255) * bug2.b) * 0.9
    bug1_power += GvB1
    print(round(bug1_power, 2))

    # Blue vs. Red
    BvR1 = (((bug1.b - ((bug1.r + bug1.g) / 10)) / 255) * bug2.r)
    bug1_power += BvR1
    print(round(bug1_power, 2))

    bug2_power = 0
    # Red vs. Green (Red is always 10% stronger)
    RvG2 = (((bug2.r - ((bug2.g + bug2.b) / 10)) / 255) * bug1.g) * 1.1
    bug2_power += RvG2
    print(round(bug2_power, 2))

    # Green vs. Blue (Blue is always 10% more defensive)
    GvB2 = (((bug2.g - ((bug2.b + bug2.r) / 10)) / 255) * bug1.b) * 0.9
    bug2_power += GvB2
    print(round(bug2_power, 2))

    # Blue vs. Red
    BvR2 = (((bug2.b - ((bug2.r + bug2.g) / 10)) / 255) * bug1.r)
    bug2_power += BvR2
    print(round(bug2_power, 2))

    return bug1_power - bug2_power


def main():
    bug1 = Bug((randint(0, 255), randint(0, 255), randint(0, 255)))
    bug2 = Bug((randint(0, 255), randint(0, 255), randint(0, 255)))



    print(bug1.rgb, bug2.rgb)
    print(fight(bug1, bug2))


if __name__ == '__main__':
    main()
