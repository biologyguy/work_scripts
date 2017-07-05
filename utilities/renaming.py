#!/usr/bin/env python3
import sys
import os
import argparse
import re
import shutil

# set up arguments for the program
parser = argparse.ArgumentParser(prog="renaming",
                                 description="A tool to help rename or delete batches of files based on RegEx queries."
                                             " No changes will happen until the -c flag is passed.")

parser.add_argument('match_string', nargs="?", help='RegEx search term')
parser.add_argument('replace_string', help='String that the RegEx will be replaced with (default is blank).', nargs="?",
                    default="")
parser.add_argument('-d', '--dir', help='Specify the directory with your files. Default is current working dir.',
                    action='store', default=os.getcwd())
parser.add_argument('-f', '--folders', help='Rename folders as well.', action='store_true', default=False)
parser.add_argument('-r', '--recursive', help='Also rename files in sub-directories.',
                    action='store_true', default=False)
parser.add_argument('-n', '--num', help='Add sequential zero-padded numbers to the end of the replacement string. '
                                        'The Int provided is the size of the padding.',
                    type=int, default=False, metavar='[INT]')
parser.add_argument('-q', '--quiet', help='Suppress all stdout output.', action='store_true', default=False)
parser.add_argument('-hd', '--hidden', help='Include hidden files in search.', action='store_true', default=False)
parser.add_argument('-o', '--overwrite', help='Allow files to be overwritten. Be careful...',
                    action='store_true', default=False)
parser.add_argument('-c', '--commit', help='Commit changes. Caution, this cannot be undone.',
                    action='store_true', default=False)
parser.add_argument('-a', '--append', help='Leave the search string intact, and append the replace string to the [f] '
                                           'front or [b] back.', type=str, choices=['f', 'b'])
parser.add_argument('-uc', '--uppercase',  action="store", type=int, choices=[1, 2, 3],
                    help='Convert first letter of first word only (1), first letter of all words (2), '
                         'or all letters (3) to upper case.')
parser.add_argument('-ns', '--num_subs', help='Only substitute regex match in file name [int] times'
                                              ' (from left to right).', type=int, default=False, metavar='[INT]')
parser.add_argument('-rm', '--remove', help='Delete files that contain a regex match',
                    action='store_true', default=False)

incoming_args = parser.parse_args()

if not incoming_args.match_string and not incoming_args.uppercase:
    sys.exit("Error: renaming requires a match string or the --uppercase flag to work. Use renaming.py -h for details.")


def stdout(string):
    if not incoming_args.quiet:
        sys.stdout.write("%s\n" % string)


def delete(file_path):
    if os.path.isdir(file_path):
        if incoming_args.commit:
            shutil.rmtree(file_path)
        output = "Delete directory (and all content)"
    else:
        if incoming_args.commit:
            os.remove(file_path)
        output = "Delete file"
    return output


def loop_names(path_list, _count):
    counter_replace = incoming_args.replace_string
    for i in path_list:
        if incoming_args.num:
            counter_replace = "%s%s" % (incoming_args.replace_string, str(_count).zfill(incoming_args.num))

        output = rename(incoming_args.match_string, i, counter_replace, dirpath)

        # output the changed files unless -q flag is passed
        if output:
            stdout("%s\t--->\t%s" % (i, output))
            _count += 1
    return _count


def rename(query, ifile, replace, path):
    # Don't include hidden files unless -hd flag is passed, and never alter . or ..
    if ifile == "." or ifile == ".." or (ifile[0] == "." and not incoming_args.hidden):
        return None

    # First deal with uppercase flag, if passed in
    if incoming_args.uppercase:
        def uc(match):
            try:
                return "%s%s" % (match.group(1), match.group(2).upper())
            except IndexError:
                return match.group(1).upper()
        output = str(ifile)
        if incoming_args.uppercase <= 2:
            output = re.sub("([a-zA-Z])", uc, output, 1)
        if incoming_args.uppercase == 2:
            output = re.sub("([^a-zA-Z.])([a-zA-Z])", uc, output)
        if incoming_args.uppercase == 3:
            output = output.upper()

    # If no fancy flags are passed, do a simple global substitution with re.sub
    elif not incoming_args.append and not incoming_args.num_subs:
        output = re.sub(query, replace, ifile)

    # Otherwise, create an iterator (i.e., when -ns or -a are passed in)
    else:
        if not any(re.finditer(query, ifile)):
            return None

        # set the number of substitutions to make, if -ns flag is passed
        num_subs = -1 if not incoming_args.num_subs else incoming_args.num_subs
        regex_iter = re.finditer(query, ifile)
        output = ifile
        offset = 0
        for i in regex_iter:
            if num_subs == 0:
                break
            num_subs -= 1

            if not incoming_args.append:
                output = "%s%s%s" % (output[:(i.start() + offset)], replace, output[(i.end() + offset):])
                if len(replace) != len(i.group()):
                    offset += len(replace) - len(i.group())
            else:
                if incoming_args.append == "f":
                    output = "%s%s%s" % (output[:(i.start() + offset)], replace, output[(i.start() + offset):])

                else:
                    output = "%s%s%s" % (output[:(i.end() + offset)], replace, output[(i.end() + offset):])

                offset += len(replace)

    if output != ifile:
        file_path = ("%s/%s" % (path, ifile))

        # If -rm flag is passed, unlink the files
        if incoming_args.remove:
            return delete(file_path)

        # Need to output a message if files will be overwritten
        if (os.path.exists("%s/%s" % (path, output)) or (output in output_list)) and not incoming_args.overwrite:
            stdout("%s not changed because %s exists (-o to overwrite)." % (ifile, output))
            return None
        if incoming_args.commit:
            os.rename(file_path, "%s/%s" % (path, output))

        output_list.append(output)
        return output
    else:
        return None

stdout("Directory: %s" % incoming_args.dir)
output_list = []
count = 1
for (dirpath, dirnames, filenames) in os.walk(incoming_args.dir):
    # always change file names
    count = loop_names(filenames, count)

    # only change directory names is -f flag is passed
    if incoming_args.folders:
        count = loop_names(dirnames, count)

    # stop in top level --dir if -r flag is not set
    if not incoming_args.recursive:
        break

if not incoming_args.commit:
    stdout("\nThese changes did not actually occur. To commit them, use the -c flag.")
