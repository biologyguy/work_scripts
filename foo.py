#!/usr/bin/env python3

import argparse
import steves_tools as br
import re
import sys
from math import floor


class Bug(object):
    def __init__(self, rgb):
        self.rgb = rgb
        self.energy = 100


def fight(bug1, bug2):
    bug1_power = 0
    # Red vs. Green (Red is always 10% stronger)
    bug1_power += (((bug1.rgb[0] - ((bug1.rgb[1] + bug1.rgb[2]) / 10)) / 255) * bug2.rgb[1]) * 1.1

    # Green vs. Blue (Blue is always 10% more defensive)
    bug1_power += (((bug1.rgb[1] - ((bug1.rgb[2] + bug1.rgb[0]) / 10)) / 255) * bug2.rgb[2]) * 0.9

    # Blue vs. Red
    bug1_power += (((bug1.rgb[2] - ((bug1.rgb[0] + bug1.rgb[1]) / 10)) / 255) * bug2.rgb[0])

    bug2_power = 0
    # Red vs. Green (Red is always 10% stronger)
    bug1_power += (((bug2.rgb[0] - ((bug2.rgb[1] + bug2.rgb[2]) / 10)) / 255) * bug1.rgb[1]) * 1.1

    # Green vs. Blue (Blue is always 10% more defensive)
    bug1_power += (((bug2.rgb[1] - ((bug2.rgb[2] + bug2.rgb[0]) / 10)) / 255) * bug1.rgb[2]) * 0.9

    # Blue vs. Red
    bug1_power += (((bug2.rgb[2] - ((bug2.rgb[0] + bug2.rgb[1]) / 10)) / 255) * bug1.rgb[0])

    return bug1_power - bug2_power


def main():
    bug1 = Bug((255, 0, 0))
    bug2 = Bug((100, 255, 0))

    print(fight(bug1, bug2))


if __name__ == '__main__':
    main()
