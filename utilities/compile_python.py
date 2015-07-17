#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Mar 19 2015 

"""
Convert a python module and dependencies into an executable.
"""

import os
import shutil
import argparse
import MyFuncs
from subprocess import Popen

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="compile_python", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("main_file", action="store",
                        help="Specify the main script to be executed.")
    parser.add_argument("-d", "--dependencies", nargs="+",
                        help="List of directories and files to be bundled with script.")
    parser.add_argument("-c", "--config",
                        help="List of dependencies in a config file. Use relative path from location of config file.")
    in_args = parser.parse_args()

    main_file = os.path.abspath(in_args.main_file)
    current_dir = os.getcwd()

    dependencies = []

    if in_args.dependencies:
        dependencies += [os.path.abspath(dpndy) for dpndy in in_args.dependencies]

    if in_args.config:
        in_args.config = os.path.abspath(in_args.config)
        os.chdir("%s" % "/".join(in_args.config.split("/")[:-1]))
        with open(in_args.config, "r") as ifile:
            configs = ifile.read().strip().split("\n")
            dependencies += [os.path.abspath(dpndy) for dpndy in configs]

    temp_dir = MyFuncs.TempDir()
    os.mkdir("%s/temp" % temp_dir.path)
    os.chdir("%s/temp" % temp_dir.path)

    script_name = str(main_file.split("/")[-1])
    script_name = ".".join(script_name.split(".")[:-1])
    shutil.copyfile(main_file,  "./__main__.py")

    if len(dependencies) > 0:
        for dpndy in dependencies:
            if os.path.isfile(dpndy):
                shutil.copyfile(dpndy, "./%s" % dpndy.split("/")[-1])
            else:
                shutil.copytree(dpndy,  "./%s" % dpndy.split("/")[-1])

    Popen("zip -rq ../new_app.zip *", shell=True).wait()
    os.chdir("../")
    Popen("echo '#!/usr/bin/env python3' | cat - new_app.zip > %s" % script_name, shell=True).wait()
    Popen("chmod +x %s" % script_name, shell=True).wait()
    os.chdir(current_dir)
    Popen("mv %s/%s ./" % (temp_dir.path, script_name), shell=True).wait()

    print('New executable created: %s' % script_name)