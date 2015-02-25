#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Feb 24 2015

"""
History recorder pushes terminal commands to a history file. For this to work, the following line must be present in
.profile or .bashrc:

    PROMPT_COMMAND="hist \$(history 1)"

Also, add a sim link to history_recorder.py called *hist* somewhere in PATH.

    $: ln -s /path/to/history_recorder.py /usr/local/bin/hist
"""



import sys
import time
from re import search, sub

date = time.strftime("%b/%d/%y %H:%M")

command = " ".join(sys.argv[2:])

# The .history lines will run if copy-pasted, but remove the datetime from the command before sending back to .history
command = sub(": [JFMASOND][a-z]{2}/[0-3][0-9]/[0-9]{2} [0-9]{2}:[0-9]{2}; ", "", command)

# I don't want .history filling up with useless crap, so skip these common 'look' commands
skip_commands = ["^ls .*", "^hs .*", "^cat .*", "^head .*", "^tail .*"]

for skip in skip_commands:
    if search(skip, command):
        sys.exit()

output = ": %s; %s" % (date, command)

current_file = "bash_history-%s" % time.strftime("%Y%m")

with open("/Volumes/Zippy/.history/%s" % current_file, "a") as ofile:
    ofile.write(bytes("%s\n" % output))