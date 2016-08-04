#!/usr/bin/env python3

# Description: access v2.topcons.net via WSDL service
# Copyright Nanjiang Shu (nanjiang.shu@scilifelab.se)

import os
import sys
import SeqBuddy as Sb
import AlignBuddy as Alb
import buddy_resources as br
import zipfile
import re
import urllib.request

try:
    from suds.client import Client
except ImportError:
    print("Please install the 'suds' package to run topcons2_wsdl.py:\n\n$ pip install suds-py3", file=sys.stderr)
    sys.exit(1)


def main(in_args):
    wsdl_url = "http://v2.topcons.net/pred/api_submitseq/?wsdl"
    fixtop = ""

    if os.path.isfile(in_args.input):
        try:
            seqbuddy = Sb.SeqBuddy(in_args.input, out_format="fasta")
            Sb.clean_seq(seqbuddy)
            # Sb.hash_ids(seqbuddy)

        except br.GuessError:
            print("Unable to read the provided input file, is it a properly formatted sequence file?")
            return 1

        if len(str(seqbuddy)) >= MAX_FILESIZE:
            print("You input seqfile is too large! Please split the file into chunks less than %d Mb."
                  % MAX_FILESIZE_IN_MB, file=sys.stderr)
            return 1

        # ***** Here's the meat ***** #
        myclient = Client(wsdl_url, cache=None)
        ret_value = myclient.service.submitjob(str(seqbuddy), fixtop, in_args.jobname, in_args.email)
        if len(ret_value) >= 1:
            jobid, result_url, numseq_str, errinfo, warninfo = ret_value[0][:5]
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
        myclient = Client(wsdl_url, cache=None)
        ret_value = myclient.service.checkjob(in_args.input)
        if len(ret_value) >= 1:
            status, result_url, errinfo = ret_value[0][:3]

            if status == "Failed":
                print("Your job with jobid %s is failed!" % in_args.input)
                if errinfo != "" and errinfo != "None":
                    print("Error message:\n%s" % errinfo)
            elif status == "Finished":
                print("Your job with jobid %s is finished!" % in_args.input)
                if not os.path.exists(in_args.outpath):
                    try:
                        os.makedirs(in_args.outpath)
                    except OSError:
                        print("Failed to create the outpath %s" % in_args.outpath)
                        return 1
                outfile = "%s/%s.zip" % (in_args.outpath, in_args.input)
                if not os.path.exists(outfile):
                    print("Retrieving")
                    urllib.request.urlretrieve(result_url, outfile)
                    if os.path.exists(outfile):
                        print("The result file %s has been retrieved for jobid %s" % (outfile, in_args.input))
                    else:
                        print("Failed to retrieve result for jobid %s" % in_args.input)

                with zipfile.ZipFile(outfile) as zf:
                    zf.extractall(in_args.outpath)

                with open("%s/%s/query.result.txt" % (in_args.outpath, in_args.input), "r") as ifile:
                    topcons = ifile.read()

                topcons = topcons.split("##############################################################################")[2:-1]

                records = []
                for rec in topcons:
                    seq_id = re.search("Sequence name: (.*)", rec).group(1).strip()
                    seq = re.search("Sequence:\n([A-Z]+)", rec).group(1).strip()
                    alignment = ""
                    for algorithm in ["TOPCONS", "OCTOPUS", "Philius", "PolyPhobius", "SCAMPI", "SPOCTOPUS"]:
                        try:
                            top_file = re.search("%s predicted topology:\n([ioMS]+)" % algorithm, rec).group(1).strip()
                            alignment += ">%s\n%s\n\n" % (algorithm, top_file)
                        except:
                            print("%s: %s" % (seq_id, algorithm))
                            pass

                    alignment = Alb.AlignBuddy(alignment)
                    # print(alignment)
                    Alb.consensus_sequence(alignment)
                    cons_seq = Sb.SeqBuddy(">%s\n%s\n" % (seq_id, seq), out_format="genbank")
                    counter = 1
                    for tmd in re.finditer("([MX]+)", str(alignment.records()[0].seq)):
                        Sb.annotate(cons_seq, "TMD%s" % counter, "%s-%s" % (tmd.start(), tmd.end()))
                        counter += 1
                    records.append(cons_seq.records[0])

                seqbuddy = Sb.SeqBuddy(records, out_format="genbank")
                seqbuddy.write("%s/%s.gb" % (in_args.outpath, in_args.input))

            elif status == "None":
                print("Your job with jobid %s does not exist! Please check you typing!" % in_args.input)
            else:
                print("Your job with jobid %s is not ready, status = %s" % (in_args.input, status))
        else:
            print("Failed to get job!")
            return 1

    return 0


MAX_FILESIZE_IN_MB = 9
MAX_FILESIZE = MAX_FILESIZE_IN_MB * 1024 * 1024

if __name__ == '__main__':
    progname = os.path.basename(sys.argv[0])
    wspace = ''.join([" "] * len(progname))

    import argparse

    parser = argparse.ArgumentParser(prog="topcons2_wsdl.py", description="Access topcons2 through a WSDL/SOAP web "
                                                                          "service (http://v2.topcons.net)")

    parser.add_argument("input", help="Location of protein input FASTA file or the job id", action="store")
    parser.add_argument("-o", "--outpath", default=os.getcwd(), action="store",
                        help="Save the retrieved data to outpath")
    parser.add_argument("-e", "--email", action="store", default="",
                        help="Specify an email address if you want a notification when the job is complete")
    parser.add_argument("-n", "--jobname", action="store", default="", help="Give the job a name")

    sys.exit(main(parser.parse_args()))
