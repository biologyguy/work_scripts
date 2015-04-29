#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Oct 8 2014

"""
Create a new template python file
"""

import argparse
import os
import datetime

parser = argparse.ArgumentParser(prog="python_templates", description="Create new python files with some useful"
                                                                      " templated code pre-written")

parser.add_argument('new_file', help='Location and name of the new file you would like to create. Defaults to current'
                                     ' working dir. All options are turned ON by default, and you can flag them off.',
                    action='store', default="%s/new_script.py" % os.getcwd())

parser.add_argument('-ns', '--no_shebang', help="Do not include location to python binary", action='store_true')
parser.add_argument('-nd', '--no_docstrings', help="Do not include docstrings", action='store_true')
parser.add_argument('-na', '--no_argparse', help='Do not include argparse module', action='store_true')
parser.add_argument('-nc', '--no_class', help="Do not include a Class block", action='store_true')
parser.add_argument('-nf', '--no_functions', help="Do not include def blocks block.", action='store_true')
parser.add_argument('-nm', '--no_main', help="Do not include a __name__ = '__main__' block.", action='store_true')
parser.add_argument('-ni', '--no_import', help="Do not import common packages", action='store_true')

parser.add_argument('-t', '--time', help="Set up a runtime counter", action='store_true')
parser.add_argument('-b', '--binary_loc', help='Specify the location of your python executable for the hashbang',
                    action='store', default="/usr/bin/env python3")
parser.add_argument('-e', '--encoding', help='Specify the encoding for the file.', action='store', default="utf-8")

in_args = parser.parse_args()


def argparse_block(white_space=0):
    out = ""
    out += "%simport argparse\n\n" % "".rjust(white_space, " ")
    out += "%sparser = argparse.ArgumentParser(prog=\"%s\", description=\"\", " \
           "formatter_class=argparse.ArgumentDefaultsHelpFormatter)\n\n" % \
           ("".rjust(white_space, " "), in_args.new_file.split("/")[-1].split(".")[0])
    out += "%sparser.add_argument(\"positional_arg1\", help=\"\", action=\"store\")\n" % "".rjust(white_space, " ")
    out += "%sparser.add_argument(\"-t\", \"--true\", help=\"\", action=\"store_true\", default=False)\n" \
           % "".rjust(white_space, " ")
    out += "%sparser.add_argument(\"-c\", \"--choice\", help=\"\", type=str, choices=[\"\", \"\"], default=False)\n" \
           % "".rjust(white_space, " ")
    out += "%sparser.add_argument(\"-m\", \"--multi_arg\", nargs=\"+\", help=\"\", default=[])\n\n" \
           % "".rjust(white_space, " ")
    out += "%sin_args = parser.parse_args()\n\n" % "".rjust(white_space, " ")
    return out

output = ""

if not in_args.no_shebang:
    output += "#!%s\n" % in_args.binary_loc

output += "# -*- coding: %s -*-\n" % in_args.encoding
today = datetime.date.today()
output += "# Created on: %s %s %s \n\n" % (today.strftime("%b"), today.day, today.year)

if not in_args.no_docstrings:
    output += '"""\nDESCRIPTION OF PROGRAM\n"""\n\n'

if not in_args.no_import:
    output += "import sys\nimport os\nimport re\nimport shutil\nimport MyFuncs\nimport SeqBuddy\n"

if in_args.time:
    output += "import timeit\n"

if not in_args.no_class:
    output += "\nclass NewClass():\n"
    if not in_args.no_docstrings:
        output += '    """DESCRIPTION OF CLASS"""\n'

    output += "    def __init__(self):\n"
    output += "        self.x = 1\n\n"
    output += "    def class_def(self):\n"
    output += "        self.x = 1\n"
    output += "        return self.x\n\n\n"

if not in_args.no_functions:
    for i in range(1, 3):
        output += "def def%s():\n" % i
        if not in_args.no_docstrings:
            output += '    """DESCRIPTION OF FUNC"""\n'

        output += "    x = 1\n"
        output += "    return x\n\n\n"


if not in_args.no_main:
    output += "if __name__ == '__main__':\n"

    if not in_args.no_argparse:
        output += argparse_block(4)

    if in_args.time:
        output += "    timer = timeit.timeit()"
    output += "    print('Hello')\n"
    if in_args.time:
        output += "    print('program executed in %s' % timeit.timeit() - timer)\n"

output = output.strip()

with open(in_args.new_file, "w") as ifile:
    ifile.write(output)

os.chmod(in_args.new_file, 0o755)