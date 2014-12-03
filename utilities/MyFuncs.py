#!/usr/bin/env python3

from multiprocessing import Process, cpu_count
from sys import stdout, exit, stderr
from time import clock
from math import floor
import os
from tempfile import TemporaryDirectory, NamedTemporaryFile
from shutil import copytree, rmtree


# might be nice to change this to a class that tracks length of last print, so it can be fully cleared when \r is called
def dynamic_print(output):
    stdout.write("\r%s\r%s" % (" " * 100, output),)
    stdout.flush()


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


def run_multicore_function(iterable, function, func_args=False, max_processes=0, quiet=False):
        # fun little piece of abstraction here... directly pass in a function that is going to be looped over, and
        # fork those loops onto independent processors. Any arguments the function needs must be provided as a list.
        cpus = cpu_count()
        if max_processes == 0:
            if cpus > 7:
                max_processes = cpus - 3
            elif cpus > 3:
                max_processes = cpus - 2
            elif cpus > 1:
                max_processes = cpus - 1
            else:
                max_processes = 1

        if max_processes > cpus:
            max_processes = cpus

        running_processes = 0
        child_list = []
        start_time = round(clock())
        elapsed = 0
        counter = 0
        if not quiet:
            print("Running function %s() on %s cores" % (function.__name__, max_processes))
        # fire up the multi-core!!
        if not quiet:
            dynamic_print("\tJob 0 of %s" % len(iterable))
    
        for next_iter in iterable:
            if type(iterable) is dict:
                next_iter = iterable[next_iter]
            while 1:     # Only fork a new process when there is a free processor.
                if running_processes < max_processes:
                    # Start new process
                    if not quiet:
                        dynamic_print("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

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
                        if (start_time + elapsed) < round(clock()):
                            elapsed = round(clock()) - start_time
                            if not quiet:
                                dynamic_print("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                        if running_processes < max_processes:
                            break

        # wait for remaining processes to complete --> this is the same code as the processor wait loop above
        if not quiet:
            dynamic_print("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

        while len(child_list) > 0:
            for i in range(len(child_list)):
                if child_list[i].is_alive():
                    continue
                else:
                    child_list.pop(i)
                    running_processes -= 1
                    break  # need to break out of the for-loop, because the child_list index is changed by pop
            if (start_time + elapsed) < round(clock()):
                elapsed = round(clock()) - start_time
                if not quiet:
                    dynamic_print("\t%s total jobs (%s, %s jobs remaining)" % (len(iterable), pretty_time(elapsed),
                                                                               len(child_list)))
            
        if not quiet:
            print(" --> DONE\n")
        # func_args = []  # This may be necessary because of weirdness in assignment of incoming arguments
        return        


class TempDir():
    def __init__(self):
        self.dir = next(self._make_dir())
        self.path = self.__str__()

    @staticmethod
    def _make_dir():
        tmp_dir = TemporaryDirectory()
        yield tmp_dir
        rmtree(tmp_dir.name)

    def save(self, location):
        if os.path.isdir(location):
            print("Save Error: Indicated output folder already exists in TempDir.save(%s)" % location, file=stderr)
            return False
        else:
            copytree(self.dir.name, location)
            return True

    def __str__(self):
        return self.dir.name


class TempFile():
    def __init__(self):
        self.handle = NamedTemporaryFile(mode='w+t')

    def write(self, content, mode="a"):
        if mode not in ["w", "a"]:
            print("Write Error: mode must be 'w' or 'a' in TempFile.write()", file=stderr)
            return False
        if mode == "a":
            self.handle.write(content)
        else:
            self.handle.truncate(0)
            self.handle.write(content)
        return True

    def home(self):
        self.handle.seek(0)
        return

    def read(self):
        self.handle.seek(0)
        return self.handle.read()

    def save(self, location):
        with open(location, "w") as ofile:
            ofile.write(self.read())
        return

    def __str__(self):
        return self.handle.name


class SafetyValve():  # Use this class if you're afraid of an infinit loop
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
            exit("Error: You just popped your global_reps safety valve. %s" % message)
    
    def test(self, state, message=""):  # test() keeps track of some variable 'state' to see if its value keeps changing
        if self.state == "%s" % state:
            self.state_reps -= 1
        else:
            self.state_reps = self._start_state_reps
            self.state = "%s" % state
            
        if self.state_reps == 0:
            exit("Error: You just popped your state_reps safety valve. %s" % message)


# Pulled this function off of Stack Overflow -- posted by nosklo
def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]
