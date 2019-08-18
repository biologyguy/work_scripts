#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: MyFuncs.py
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Collection of useful classes and functions
"""

from multiprocessing import Process, cpu_count
import sys
from time import time, sleep
from math import floor, ceil
import os
from tempfile import TemporaryDirectory
from shutil import copytree, rmtree, copyfile
import re
import string
from random import choice
import signal
import argparse
from collections import OrderedDict
from random import random
import io


# Credit to rr- (http://stackoverflow.com/users/2016221/rr)
# http://stackoverflow.com/questions/18275023/dont-show-long-options-twice-in-print-help-from-argparse
class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


class Timer(object):
    def __init__(self):
        self.start = round(time())
        self.split_time = round(time())

    def split(self):
        split = round(time()) - self.split_time
        self.split_time = round(time())
        return pretty_time(split)

    def total_elapsed(self):
        return pretty_time(round(time()) - self.start)


class RunTime(object):
    def __init__(self, prefix=None, postfix=None, out_type=sys.stdout, _sleep=0, final_clear=False):
        self.out_type = out_type
        self.prefix = prefix if prefix else ""
        self.postfix = postfix if postfix else ""
        self.running_process = None
        self.sleep = _sleep
        self.final_clear = final_clear

    def _run(self, check_file_path):
        d_print = DynamicPrint(self.out_type)
        start_time = round(time())
        elapsed = 0
        while True:
            prefix = self.prefix if not hasattr(self.prefix, '__call__') else self.prefix()
            postfix = self.postfix if not hasattr(self.postfix, '__call__') else self.postfix()
            with open("%s" % check_file_path, "r") as ifile:
                if round(time()) - start_time == elapsed:
                    continue
                elif ifile.read() == "Running":
                    d_print.write("%s%s%s" % (prefix, pretty_time(elapsed), postfix))
                    elapsed = round(time()) - start_time
                else:
                    if not self.final_clear:
                        d_print.write("%s%s%s\n" % (prefix, pretty_time(elapsed), postfix))
                    else:
                        d_print.clear()
                    break
            sleep(self.sleep)
        return

    def start(self):
        if self.running_process:
            self.end()
        tmp_file = TempFile()
        tmp_file.write("Running")
        p = Process(target=self._run, args=(tmp_file.path,))
        p.daemon = 1
        p.start()
        self.running_process = [tmp_file, p]
        return

    def end(self):
        if not self.running_process:
            return
        self.running_process[0].clear()
        while self.running_process[1].is_alive():
            continue
        self.running_process = None
        return


# maybe use curses library in the future to extend this for multi-line printing
class DynamicPrint(object):
    def __init__(self, out_type="stdout", quiet=False):
        self._last_print = ""
        self._next_print = ""
        self._writer = self._write()

        out_type = sys.stdout if out_type == "stdout" else out_type
        out_type = sys.stderr if out_type == "stderr" else out_type
        self.out_type = out_type
        self.quiet = quiet

    def _write(self):
        try:
            while True:
                self.out_type.write("\r%s\r%s" % (" " * len(self._last_print), self._next_print),)
                self.out_type.flush()
                self._last_print = self._next_print
                yield
        finally:
            self.out_type.write("")

    def write(self, content):
        if not self.quiet:
            content = re.sub("\t", "    ", content)
            self._next_print = content
            next(self._writer)
        return

    def new_line(self, number=1):
        if not self.quiet:
            self.out_type.write("\n" * number)
            self.out_type.flush()
            self._last_print = ""
        return

    def clear(self):
        self.write("")
        return


class KellysColors(object):
    # https://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
    def __init__(self, convert2hex=False):
        self.kelly_colors = OrderedDict(deep_yellowish_brown=(89, 51, 21),
                                        strong_reddish_brown=(127, 24, 13),
                                        strong_purplish_red=(179, 40, 81),
                                        strong_purplish_pink=(246, 118, 142),
                                        vivid_red=(193, 0, 32),
                                        vivid_reddish_orange=(241, 58, 19),
                                        vivid_orange=(255, 104, 0),
                                        strong_yellowish_pink=(255, 122, 92),
                                        vivid_orange_yellow=(255, 142, 0),
                                        vivid_yellow=(255, 179, 0),
                                        vivid_greenish_yellow=(244, 200, 0),
                                        grayish_yellow=(206, 162, 98),
                                        vivid_yellowish_green=(147, 170, 0),
                                        vivid_green=(0, 125, 52),
                                        dark_olive_green=(35, 44, 22),
                                        very_light_blue=(166, 189, 215),
                                        strong_blue=(0, 83, 138),
                                        strong_violet=(83, 55, 122),
                                        strong_purple=(128, 62, 117),
                                        medium_gray=(129, 112, 102))
        self.hex = convert2hex
        self.generator = self.color_iter()

    def color_iter(self):
        degree = 0
        while True:
            for color, rgb in self.kelly_colors.items():
                rgb = purturb_rgb(rgb, degree)
                if self.hex:
                    rgb = '#%02x%02x%02x' % rgb
                yield rgb
            degree += 10

    def iter(self):
        return next(self.generator)


def purturb_rgb(rgb, degree=10):
    new_rgb = []
    for code in rgb:
        degree = round(random() * degree)
        degree = degree * -1 if random() < 0.5 else degree
        while degree != 0:
            code += degree
            if code > 255:
                degree = 255 - code
                code = 255
            elif code < 0:
                degree = abs(code)
                code = 0
            else:
                degree = 0
        new_rgb.append(code)
    return new_rgb[0], new_rgb[1], new_rgb[2]


def pretty_time(seconds):
    if seconds < 60:
        output = "%i sec" % seconds
    elif seconds < 3600:
        minutes = floor(seconds / 60)
        seconds -= minutes * 60
        output = "%i min, %i sec" % (minutes, seconds)
    elif seconds < 86400:
        hours = floor((seconds / 60) / 60)
        seconds -= hours * 60 * 60
        minutes = floor(seconds / 60)
        seconds -= minutes * 60
        output = "%i hrs, %i min, %i sec" % (hours, minutes, seconds)
    else:
        days = floor(((seconds / 60) / 60) / 24)
        seconds -= (days * 60 * 60 * 24)
        hours = floor((seconds / 60) / 60)
        seconds -= (hours * 60 * 60)
        minutes = floor(seconds / 60)
        seconds -= (minutes * 60)
        output = "%i days, %i hrs, %i min, %i sec" % (days, hours, minutes, seconds)

    return output


def pretty_number(num, mode='short', precision=2):  # mode in ['short', 'medium', 'long']
    magnitude = 0
    while abs(num) >= 1000 and magnitude < 8:
        magnitude += 1
        num = round(num / 1000.0, precision)
    if mode == 'short':
        return '%s %s' % (num, ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'][magnitude])
    elif mode == 'medium':
        return '%s %s' % (num, ['', 'Kilo', 'Mega', 'Giga', 'Tera', 'Peta', 'Exa', 'Zetta', 'Yotta'][magnitude])
    elif mode == 'long':
        return '%s %s' % (num, ['', 'Thousand', 'Million', 'Billion', 'Trillion', 'Quadrillion', 'Quintillion',
                                'Sextillion', 'Septillion'][magnitude])
    else:
        raise ValueError("Valid 'mode' values are 'short', 'medium', and 'long'")


def usable_cpu_count():
    cpus = cpu_count()
    if cpus > 7:
        max_processes = cpus - 3
    elif cpus > 3:
        max_processes = cpus - 2
    elif cpus > 1:
        max_processes = cpus - 1
    else:
        max_processes = 1

    return max_processes


def run_multicore_function(iterable, function, func_args=False, max_processes=0, quiet=False, out_type=sys.stdout):
    """
    fun little piece of abstraction here... directly pass in a function that is going to be looped over, and
    fork those loops onto independent processes. Any arguments the function needs must be provided as a list.
    :param iterable: Iterable collection of things that will be fed into the multicore function
    :param function: Reference to a function that each iterable item will be fed into
    :param func_args: Tuple of any additional arguments that need to be passed into the function
    :param max_processes: Either auto-detect the system resources (set 0) or specify a number
    :param quiet: Do not print status to console
    :param out_type: Switch output to stderr, if desired.
    :return:
    """
    d_print = DynamicPrint(out_type)

    if max_processes == 0:
        max_processes = usable_cpu_count()

    else:
        cpus = cpu_count()
        if max_processes > cpus:
            max_processes = cpus
        elif max_processes < 1:
            max_processes = 1

    max_processes = max_processes if max_processes < len(iterable) else len(iterable)

    running_processes = 0
    child_list = []
    start_time = round(time())
    elapsed = 0
    counter = 0
    if not quiet:
        d_print.write("Running function %s() on %s cores\n" % (function.__name__, max_processes))
    # fire up the multi-core!!
    if not quiet:
        d_print.write("\tJob 0 of %s" % len(iterable))

    for next_iter in iterable:
        if type(iterable) is dict:
            next_iter = iterable[next_iter]
        while 1:     # Only fork a new process when there is a free processor.
            if running_processes < max_processes:
                # Start new process
                if not quiet:
                    d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                if func_args:
                    if not isinstance(func_args, list):
                        exit("Error in run_multicore_function(): The arguments passed into the multi-thread "
                             "function must be provided as a list")
                    p = Process(target=function, args=(next_iter, func_args))

                else:
                    p = Process(target=function, args=(next_iter,))
                p.start()
                child_list.append(p)
                running_processes += 1
                counter += 1
                break
            else:
                # processor wait loop
                while 1:
                    for i in range(len(child_list)):
                        if child_list[i].is_alive():
                            continue
                        else:
                            child_list.pop(i)
                            running_processes -= 1
                            break

                    if not quiet:
                        if (start_time + elapsed) < round(time()):
                            elapsed = round(time()) - start_time
                            d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                    if running_processes < max_processes:
                        break

    # wait for remaining processes to complete --> this is the same code as the processor wait loop above
    if not quiet:
        d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

    while len(child_list) > 0:
        for i in range(len(child_list)):
            if child_list[i].is_alive():
                continue
            else:
                child_list.pop(i)
                running_processes -= 1
                break  # need to break out of the for-loop, because the child_list index is changed by pop

        if not quiet:
            if (start_time + elapsed) < round(time()):
                elapsed = round(time()) - start_time
                d_print.write("\t%s total jobs (%s, %s jobs remaining)" % (len(iterable), pretty_time(elapsed),
                                                                           len(child_list)))

    if not quiet:
        d_print.write("\tDONE: %s jobs in %s\n" % (len(iterable), pretty_time(elapsed)))
    # func_args = []  # This may be necessary because of weirdness in assignment of incoming arguments
    return


class TempDir(object):
    def __init__(self):
        self.dir = next(self._make_dir())
        self.path = self.dir.name
        self.subdirs = []
        self.subfiles = []

    def _make_dir(self):
        tmp_dir = TemporaryDirectory()
        yield tmp_dir
        rmtree(self.path)

    def copy_to(self, src):
        full_path = os.path.abspath(src)
        end_path = os.path.split(full_path)[1]
        if os.path.isdir(src):
            copytree(src, os.path.join(self.path, end_path))
        elif os.path.isfile(src):
            copyfile(src, os.path.join(self.path, end_path))
        else:
            return False
        return os.path.join(self.path, end_path)

    def subdir(self, dir_name=None):
        if not dir_name:
            dir_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            while dir_name in self.subdirs:  # Catch the very unlikely case that a duplicate occurs
                dir_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])

        subdir_path = os.path.join(self.path, dir_name)
        if not os.path.exists(subdir_path):
            os.mkdir(subdir_path)
        if dir_name not in self.subdirs:
            self.subdirs.append(dir_name)
        return subdir_path

    def del_subdir(self, _dir):
        path, _dir = os.path.split(_dir)
        del self.subdirs[self.subdirs.index(_dir)]
        rmtree(os.path.join(self.path, _dir))
        return

    def subfile(self, file_name=None):
        if not file_name:
            file_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])
            root, dirs, files = next(walklevel(self.path))
            while file_name in files:  # Catch the very unlikely case that a duplicate occurs
                file_name = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])

        open(os.path.join(self.path, file_name), "w", encoding="utf-8").close()
        self.subfiles.append(file_name)
        return os.path.join(self.path, file_name)

    def del_subfile(self, _file):
        path, _file = os.path.split(_file)
        del self.subfiles[self.subfiles.index(_file)]
        os.remove(os.path.join(self.path, _file))
        return

    def save(self, location, keep_hash=False):
        location = location if not keep_hash else os.path.join(location, os.path.split(self.path)[-1])
        if os.path.isdir(location):
            print("Save Error: Indicated output folder already exists in TempDir.save(%s)" % location, file=sys.stderr)
            return False
        else:
            copytree(self.dir.name, location)
            return True


class TempFile(object):
    # I really don't like the behavior of tempfile.[Named]TemporaryFile(), so hack TemporaryDirectory() via TempDir()
    def __init__(self, mode="w", byte_mode=False, encoding="utf-8"):
        self._tmp_dir = TempDir()  # This needs to be a persistent (ie self.) variable, or the directory will be deleted
        path, dir_hash = os.path.split(self._tmp_dir.path)
        self.name = dir_hash
        self.path = os.path.join(self._tmp_dir.path, dir_hash)
        open(self.path, "w", encoding=encoding).close()
        self.handle = None
        self.bm = "b" if byte_mode else ""
        self.mode = mode
        self.encoding = encoding

    def open(self, mode=None):
        mode = "%s%s" % (self.mode, self.bm) if not mode else "%s%s" % (mode, self.bm)
        if self.handle:
            self.close()
        encoding = None if self.bm or "b" in mode else self.encoding
        self.handle = open(self.path, mode, encoding=encoding)

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None

    def get_handle(self, mode=None):
        self.open(self.mode) if not mode else self.open(mode)
        return self.handle

    def write(self, content, mode="a"):
        mode = "%s%s" % (mode, self.bm)
        if mode not in ["w", "wb", "a", "ab"]:
            print("Write Error: mode must be 'w' or 'a' in TempFile.write()", file=sys.stderr)
            return False
        already_open = True if self.handle else False
        if not already_open:
            self.open(mode[0])
        if mode in ["a", "ab"]:
            self.handle.write(content)
        else:
            self.handle.truncate(0)
            self.handle.write(content)
        if not already_open:
            self.close()
        return True

    def read(self):
        already_open = True if self.handle else False
        position = 0
        if already_open:
            position = self.handle.tell()
            self.close()
        encoding = None if self.bm else self.encoding
        with open(self.path, "r%s" % self.bm, encoding=encoding) as ifile:
            content = ifile.read()
        if already_open:
            self.open(mode="a")
            self.handle.seek(position)
        return content

    def clear(self):
        self.close()
        content = "" if self.bm == "" else b""
        self.write(content, mode="w")
        return

    def save(self, location):
        encoding = None if self.bm else self.encoding
        with open(location, "w%s" % self.bm, encoding=encoding) as ofile:
            ofile.write(self.read())
        return


class SafetyValve(object):  # Use this class if you're afraid of an infinite loop
    def __init__(self, global_reps=1000, state_reps=10, counter=0):
        self.counter = counter

        self._start_global_reps = global_reps
        self.global_reps = global_reps

        self._start_state_reps = state_reps
        self.state_reps = state_reps

        self.state = ""

    def step(self, message=""):  # step() is general, and doesn't care whether useful computation is on going
        self.global_reps -= 1
        self.counter += 1
        if self.global_reps == 0:
            raise RuntimeError("You just popped your global_reps safety valve. %s" % message)
        return True

    def test(self, state, message=""):  # test() keeps track of some variable 'state' to see if its value keeps changing
        if self.state == str(state):
            self.state_reps -= 1
        else:
            self.state_reps = self._start_state_reps
            self.state = str(state)

        if self.state_reps == 0:
            raise RuntimeError("You just popped your state_reps safety valve. %s" % message)
        return True


# Pulled this function off of Stack Overflow -- posted by nosklo
# Note that this is a generator, so need to use next() or `with` to get a result
def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]


def copydir(source, dest):
    for root, dirs, files in os.walk(source):
        if not os.path.isdir(dest):
            os.makedirs(dest)
        for each_file in files:
            rel_path = re.sub(source, '', root)
            rel_path = rel_path.lstrip(os.sep)
            dest_path = os.path.join(dest, rel_path, each_file)
            copyfile(os.path.join(root, each_file), dest_path)


def normalize(data, trim_ends=1.0):
    if 0. > trim_ends > 1.0:
        raise ValueError("normalize() trim_ends parameter should be between 0.5 and 1.0")

    if trim_ends > 0.5:
        trim_ends = 1 - trim_ends

    max_limit = ceil(len(data) * (1 - trim_ends)) - 1
    min_limit = -1 * (max_limit + 1)

    sorted_data = sorted([data[key] for key in data]) if type(data) == dict else sorted(data)
    _max = sorted_data[max_limit]
    _min = sorted_data[min_limit]
    data_range = _max - _min

    keys = [key for key in data] if type(data) == dict else range(len(data))
    for key in keys:
        data[key] = (data[key] - _min) / data_range
        data[key] = 1. if data[key] > 1. else data[key]
        data[key] = 0. if data[key] < 0. else data[key]

    return data


# This will only work if SMTP is running locally
def sendmail(sender, recipient, subject, message):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    msg = MIMEMultipart()
    msg.preamble = subject
    msg.add_header("From", sender)
    msg.add_header("Subject", subject)
    msg.add_header("To", recipient)

    msg.attach(MIMEText(message))

    smtp = smtplib.SMTP('localhost')
    smtp.starttls()
    smtp.sendmail(sender, recipient, msg.as_string())
    smtp.quit()
    return


def ask(input_prompt, default="yes", timeout=0):
    if default == "yes":
        yes_list = ["yes", "y", '']
        no_list = ["no", "n", "abort"]
    else:
        yes_list = ["yes", "y"]
        no_list = ["no", "n", "abort", '']

    def kill(*args):
        raise TimeoutError

    try:
        signal.signal(signal.SIGALRM, kill)
        signal.alarm(timeout)
        _response = input(input_prompt)
        signal.alarm(0)
        while True:
            if _response.lower() in yes_list:
                return True
            elif _response.lower() in no_list:
                return False
            else:
                print("Response not understood. Valid options are 'yes' and 'no'.")
                signal.alarm(timeout)
                _response = input(input_prompt)
                signal.alarm(0)

    except TimeoutError:
        return False


def tablizer(raw_input, max_width=110, cell_widths=(), right_pad=2, underline_heading=True, to_list=False):
    """
    Convert a list of lists of strings into a nice table
    :param raw_input: [(row1, row2, row3, ...), (row1, row2, row3, ...)]
    :param max_width: Maximum width of the entire table
    :param cell_widths: Granular control of how big each cell can be
    :param right_pad: How much space between cells?
    :param underline_heading: Add underlines to table headers
    :param to_list: Output the result as a list, instead of a string
    :return: str
    """
    # Start by getting the incoming input into a list of lists of strings
    raw_input = [[str(item) for item in group] for group in raw_input]

    if not cell_widths:
        max_widths = [[indx, 0] for indx in range(len(raw_input[0]))]
        largest_words = [[indx, ""] for indx in range(len(raw_input[0]))]
        for row in raw_input:
            for indx, item in enumerate(row):
                item = str(item)
                max_widths[indx][1] = len(item) if len(item) > max_widths[indx][1] else max_widths[indx][1]
                if not item or re.match(" +$", item):
                    words_sizes = [" "]
                else:
                    words_sizes = sorted(item.split(), key=lambda x: len(x), reverse=True)
                largest_words[indx][1] = words_sizes[0] if len(words_sizes[0]) > len(largest_words[indx][1]) \
                    else largest_words[indx][1]

        # This is a check that max_width is sufficient to cover the largest words in each column. Expand if needed
        min_width = sum([len(x[1]) for x in largest_words]) + (right_pad * len(raw_input[0]))
        max_width = max_width if max_width > min_width else min_width

        max_widths = sorted(max_widths, key=lambda x: x[1])
        remaining_cols = len(raw_input[0])
        remaining_space = max_width - (remaining_cols * right_pad)
        equal_width = floor(remaining_space / remaining_cols)
        cell_widths = [0 for _ in range(remaining_cols)]
        for indx, val in max_widths:
            if val + right_pad < equal_width:
                cell_widths[indx] = val + right_pad
                remaining_cols -= 1
                remaining_space -= (val + right_pad)
                equal_width = remaining_space if not remaining_cols else floor(remaining_space / remaining_cols)
            else:
                cell_widths[indx] = equal_width + right_pad
    for col_indx, mw in enumerate(cell_widths):
        for raw_row in raw_input:
            raw_row[col_indx] = chunk_text(raw_row[col_indx], mw).split("\n")

    output = "" if not to_list else []

    for row_indx, row in enumerate(raw_input):
        num_lines = max([len(x) for x in row])
        for line_indx in range(num_lines):
            line = ""
            for cell_indx, lines in enumerate(row):
                if len(lines) >= line_indx + 1:
                    if underline_heading and row_indx == 0:
                        line += "\033[4m"
                        line += lines[line_indx]
                        line += "\033[24m"
                        line += " " * (cell_widths[cell_indx] - len(lines[line_indx]))
                    else:
                        line += lines[line_indx].ljust(cell_widths[cell_indx], " ")
                else:
                    line += " ".ljust(cell_widths[cell_indx], " ")
            output += line + "\n" if not to_list else [line + "\n"]
    return output


def chunk_list(l, num_chunks):
    """
    Break up a list into a list of lists
    :param l: Input list
    :param num_chunks: How many lists should the list be chunked into
    :return:
    """
    num_chunks = int(num_chunks)
    if num_chunks < 1 or not l:
        raise AttributeError("Input list must have items in it and num_chunks must be a positive integer")

    size = int(ceil(len(l) / num_chunks))
    num_long = len(l) % num_chunks
    num_long = num_long if num_long != 0 else num_chunks
    chunks = [l[i:i + size] for i in range(0, num_long * size, size)]
    if size != 1:
        chunks += [l[i:i + size - 1] for i in range(num_long * size, len(l), size - 1)]
    return chunks


def chunk_text(text, max_line_len):
    output = ""
    text = text.split()
    cur_line = ""
    for word in text:
        if len(cur_line) + len(word) + 1 > max_line_len:
            output += cur_line + "\n"
            cur_line = word
        else:
            cur_line += " " + word
    output += cur_line
    return output


def easy_read(ifile):
    if str(type(ifile)) == "<class '_io.TextIOWrapper'>":
        if not ifile.seekable():
            input_txt = ifile.read()
            ifile = io.StringIO(utf_encode(input_txt))
        ifile = ifile.read()
    elif os.path.isfile(ifile):
        in_file = ifile
        with open(in_file, "r", encoding="utf-8") as ifile:
            ifile = io.StringIO(ifile.read())
            ifile = ifile.read()
    else:
        ifile = None
    return ifile


def utf_encode(data):
    import codecs
    tmp = TempFile()
    with open(tmp.path, "w", encoding="utf-8") as ofile:
        ofile.write(data)

    with codecs.open(tmp.path, "r", "utf-8", errors="replace") as ifile:
        data = ifile.read()
    data = re.sub(r'\r', r'\n', data)
    return data


def num_sorted(input_list):
    """
    Sort a list of strings in the way that takes embedded numbers into account
    """
    def convert(text):
        return int(text) if text.isdigit() else text

    def alpha_num_key(key):
        return [convert(c) for c in re.split('([0-9]+)', str(key))]

    return sorted(input_list, key=alpha_num_key)
