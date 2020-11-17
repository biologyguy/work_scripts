#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Feb 24 2015

"""
Add a sim link to history_recorder.py called *hist* somewhere in $PATH.

    $: ln -s /path/to/history_recorder.py /usr/local/bin/hist

History recorder pushes terminal commands to a history file. For this to work, the following lines must be present in
.profile or .bashrc:

    export HISTORYREC=/path/to/home   # Or to some other location on the system
    PROMPT_COMMAND="hist $HISTORYREC \$(history 1) $$"

Make a directory called .history

    mkdir $HISTORYREC/.history
"""

import sys
import time
from re import search, sub

date = time.strftime("%b/%d/%y %H:%M")

command = " ".join(sys.argv[3:-1])
pid = sys.argv[-1]

# The .history lines will run if copy-pasted, but remove the datetime from the command before sending back to .history
command = sub(":[0-9]* [JFMASOND][a-z]{2}/[0-3][0-9]/[0-9]{2} [0-9]{2}:[0-9]{2}; ", "", command)

# I don't want .history filling up with useless crap, so skip these common 'look' commands
skip_commands = ["ls", "la", "ll", "hs", "js", "cat", "head", "tail"]

for skip in skip_commands:
    if search(r'^%s *' % skip, command) and not search(r'\|', command):
        sys.exit()

# To override skips, prepend :; to the front of command, and then strip those characters from the history
if command[:2] == ":;":
    command = command[2:].strip()

output = ":%s %s; %s" % (pid, date, command)

current_file = "bash_history-%s" % time.strftime("%Y%m")

with open("%s/.history/%s" % (sys.argv[1], current_file), "a") as ofile:
    ofile.write("%s\n" % output)
