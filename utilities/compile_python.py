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
                        help="Specify the main script to be executed. All files and directories in the same directory "
                             "as that file will be included in the final executable.")

    in_args = parser.parse_args()

    temp_dir = MyFuncs.TempDir()
    current_dir = os.getcwd()
    main_dir = os.path.abspath(in_args.main_file).split("/")
    script_name = str(main_dir[-1])
    main_dir = "/".join(main_dir[:-1])
    shutil.copytree(main_dir, "%s/temp" % temp_dir.path)
    os.chdir("%s/temp" % temp_dir.path)
    os.rename(script_name, "__main__.py")

    Popen("zip -r ../new_app.zip *", shell=True).wait()
    os.chdir(temp_dir.path)
    script_name = ".".join(script_name.split(".")[:-1])
    Popen("echo '#!/usr/bin/env python3' | cat - new_app.zip > %s" % script_name, shell=True).wait()
    Popen("chmod +x %s" % script_name, shell=True).wait()
    os.rename(script_name, "%s/%s" % (main_dir, script_name))
    os.chdir(current_dir)

    print('New executable created: %s' % script_name)