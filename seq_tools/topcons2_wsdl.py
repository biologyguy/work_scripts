#!/usr/bin/env python3

# Description: access v2.topcons.net via WSDL service
# Copyright Nanjiang Shu (nanjiang.shu@scilifelab.se)

import os
import sys
progname = os.path.basename(sys.argv[0])
wspace = ''.join([" "] * len(progname))

try:
    from suds.client import Client
except ImportError:
    print("Please install the 'suds' package to run topcons2_wsdl.py:\n\n$ pip install suds-py3", file=sys.stderr)
    sys.exit(1)

import urllib.request

MAX_FILESIZE_IN_MB = 9
MAX_FILESIZE = MAX_FILESIZE_IN_MB * 1024 * 1024

usage_short = """
Usage: %s -m submit|get [-seq SEQFILE] [-jobname NAME] [-email EMAIL]
       %s               [-fix FIXFILE]
       %s               [-jobid JOBID] [-outpath DIR]
""" % (progname, wspace, wspace)

usage_ext="""
Description:
    Access topcons2 web-server (http://v2.topcons.net) through WSDL service

OPTIONS:
  -m submit|get  Set the mode
                 submit - submit a job to WSDL
                 get    - retrieve the result from the server

  -seq    FILE   Supply input sequence in FASTA format

  -jobname STR   Give the job a name

  -email   STR   Send a notification to the email when the result is ready

  -jobid   STR   Retrieve the result by supplying a valid jobid

  -outpath DIR   Save the retrieved data to outpath, (default: ./)

  -fix    FILE   Supply a file with topology constraints (TO BE IMPLEMENTED)

  -h, --help     Print this help message and exit

Created 2015-02-04, updated 2015-02-06, Nanjiang Shu
"""

usage_exp = """
Examples:
    # submit test.fa with jobname 'test' to the server
    %s -m submit -seq test.fa -jobname test

    # try to retrieve the result for jobid 'rst_TTT' and save it to the current
    # directory
    %s -m get -jobid rst_TTT

""" % (progname, progname)


def my_getopt_str(argv, i):
    """
    Get a string from the argument list, return the string and the updated
    index to the argument list
    """
    try:
        opt = argv[i + 1]
        if opt[0] == "-":
            msg = "Error! option '%s' must be followed by a string, not an option arg."
            print(msg % argv[i], file=sys.stderr)
            sys.exit(1)
        return opt, i + 2
    except IndexError:
        msg = "Error! option '%s' must be followed by a string"
        print(msg % argv[i], file=sys.stderr)
        sys.exit(1)


def print_help(fpout=sys.stdout):
    print(usage_short, file=fpout)
    print(usage_ext, file=fpout)
    print(usage_exp, file=fpout)


def read_file(infile, mode="r"):
    try:
        fpin = open(infile, mode)
        content = fpin.read()
        fpin.close()
        return content
    except IOError:
        print("Failed to read file %s with mode '%s'" % (infile, mode), file=sys.stderr)
        return ""


def main(g_params):
    argv = sys.argv
    num_argv = len(argv)
    if num_argv < 2:
        print_help()
        return 1

    wsdl_url = "http://v2.topcons.net/pred/api_submitseq/?wsdl"
    mode = ""
    jobid = ""
    email = ""
    jobname = ""
    fixtopfile = ""
    seqfile = ""
    outpath = "./"

    i = 1
    is_non_option_arg = False
    while i < num_argv:
        if is_non_option_arg:
            print("Error! Wrong argument: %s" % argv[i], file=sys.stderr)
            return 1

        elif argv[i] == "--":
            is_non_option_arg = True
            i += 1

        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                print_help()
                return 0
            elif argv[i] in ["-m", "--m"]:
                (mode, i) = my_getopt_str(argv, i)
            elif argv[i] in ["-seq", "--seq"]:
                (seqfile, i) = my_getopt_str(argv, i)
            elif argv[i] in ["-jobname", "--jobname"]:
                (jobname, i) = my_getopt_str(argv, i)
            elif argv[i] in ["-email", "--email"]:
                (email, i) = my_getopt_str(argv, i)
            elif argv[i] in ["-fix", "--fix"]:
                (fixtopfile, i) = my_getopt_str(argv, i)
            elif argv[i] in ["-jobid", "--jobid"]:
                (jobid, i) = my_getopt_str(argv, i)
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = my_getopt_str(argv, i)
            else:
                print("Error! Wrong argument: %s" % argv[i], file=sys.stderr)
                return 1
        else:
            print("Error! Wrong argument: %s" % argv[i], file=sys.stderr)
            return 1

    if mode == "":
        print("mode not set. exit!", file=sys.stderr)
        print(usage_short)
        return 1
    elif mode not in ["submit", "get"]:
        print("unrecognized mode. exit!", file=sys.stderr)
        print(usage_short)
        return 1

    if mode == "submit":
        if seqfile == "":
            print("You want to submit a job but seqfile not set. exit!", file=sys.stderr)
            print(usage_short)
            return 1
        elif not os.path.exists(seqfile):
            print("seqfile %s does not exist. exit!" % seqfile, file=sys.stderr)
            return 1

        try:
            filesize = os.path.getsize(seqfile)
        except OSError:
            print("failed to get the size of seqfile %s. exit" % seqfile, file=sys.stderr)
            return 1

        if filesize >= MAX_FILESIZE:
            print("You input seqfile %s exceeds the upper limit %d Mb." % (seqfile, MAX_FILESIZE_IN_MB),
                  file=sys.stderr)
            print("Please split your seqfile and submit again.", file=sys.stderr)
            return 1
        seq = read_file(seqfile)

        fixtop = ""
        if fixtopfile != "":
            fixtop = read_file(fixtopfile)
        myclient = Client(wsdl_url, cache=None)
        ret_value = myclient.service.submitjob(seq, fixtop, jobname, email)
        if len(ret_value) >= 1:
            strs = ret_value[0]
            jobid = strs[0]
            # result_url = strs[1]
            numseq_str = strs[2]
            errinfo = strs[3]
            warninfo = strs[4]
            if jobid != "None" and jobid != "":
                print("You have successfully submitted your job with %s sequences. jobid = %s" % (numseq_str, jobid))
                if warninfo != "" and warninfo != "None":
                    print("Warning message:\n%s" % warninfo)
            else:
                print("Failed to submit job!")
                if errinfo != "" and errinfo != "None":
                    print("Error message:\n%s" % errinfo)
                if warninfo != "" and warninfo != "None":
                    print("Warning message:\n%s" % warninfo)
        else:
            print("Failed to submit job!")
            return 1
    else:
        if jobid == "":
            print("You want to get the result of a job but jobid not set. exit!", file=sys.stderr)
            return 1
        myclient = Client(wsdl_url, cache=None)
        ret_value = myclient.service.checkjob(jobid)
        if len(ret_value) >= 1:
            strs = ret_value[0]
            status = strs[0]
            result_url = strs[1]
            errinfo = strs[2]
            if status == "Failed":
                print("Your job with jobid %s is failed!" % jobid)
                if errinfo != "" and errinfo != "None":
                    print("Error message:\n%s" % errinfo)
            elif status == "Finished":
                print("Your job with jobid %s is finished!" % jobid)
                if not os.path.exists(outpath):
                    try:
                        os.makedirs(outpath)
                    except OSError:
                        print("Failed to create the outpath %s" % outpath)
                        return 1
                outfile = "%s/%s.zip" % (outpath, jobid)
                urllib.request.urlretrieve(result_url, outfile)
                if os.path.exists(outfile):
                    print("The result file %s has been retrieved for jobid %s" % (outfile, jobid))
                else:
                    print("Failed to retrieve result for jobid %s" % jobid)
            elif status == "None":
                print("Your job with jobid %s does not exist! Please check you typing!" % jobid)
            else:
                print("Your job with jobid %s is not ready, status = %s" % (jobid, status))
        else:
            print("Failed to get job!")
            return 1

    return 0


def init_global_parameter():
    return {'isQuiet': True}

if __name__ == '__main__':
    in_args = init_global_parameter()
    sys.exit(main(in_args))
