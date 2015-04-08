#!/usr/bin/python

from subprocess import Popen

#I was porting this from bash into python, but haven't finished yet. 

# This is a simple script which will carry out all of the basic steps
# required to make a PSIPRED V2 prediction. Note that it assumes that the
# following programs are in the appropriate directories:
# blastpgp - PSIBLAST executable (from NCBI toolkit)
# makemat - IMPALA utility (from NCBI toolkit)
# psipred - PSIPRED V3 program
# psipass2 - PSIPRED V3 program

# NOTE: Script modified to be more cluster friendly (DTJ April 2008)

#argeParse stuff
parser = argparse.ArgumentParser(prog="psipred", description="wrapper");
parser.add_argument('fasta', help='');
parser.add_argument('-d','--dbname', default='/usr/local/blastdbs/uniref90/uniref90_filtered', action='store', help='');
parser.add_argument('out_dir', help='', action='store');
in_args = parser.parse_args()

# The name of the BLAST data bank
dbname = "/usr/local/blastdbs/uniref90/uniref90_filtered"

# Where the PSIPRED V2 programs have been installed
execdir = "./bin"

# Where the PSIPRED V2 data files have been installed
datadir = "./data"

basename = in_args.fasta.split("/")[-1].split(".")[0]
rootname = $basename:t

# Generate a "unique" temporary filename root
hostid = `hostid`
tmproot = psitmp$$$hostid

\cp -f $1 $tmproot.fasta

print "Running PSI-BLAST with sequence" + in_args.fasta

Popen("blastpgp -b 0 -j 3 -h 0.001 -d " + in_args.dbname + " -i " + in_args.fasta + " -C " + $tmproot.chk + " >& " + $tmproot.blast,shell=True).wait()

if ($status != 0) then
    tail $tmproot.blast
    echo "FATAL: Error whilst running blastpgp - script terminated!"
    exit $status
endif

echo "Predicting secondary structure..."

echo $tmproot.chk > $tmproot.pn
echo $tmproot.fasta > $tmproot.sn

$ncbidir/makemat -P $tmproot

if ($status != 0) then
    echo "FATAL: Error whilst running makemat - script terminated!"
    exit $status
endif

echo Pass1 ...

$execdir/psipred $tmproot.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $rootname.ss

if ($status != 0) then
    echo "FATAL: Error whilst running psipred - script terminated!"
    exit $status
endif

echo Pass2 ...

$execdir/psipass2 $datadir/weights_p2.dat 1 1.0 1.0 $rootname.ss2 $rootname.ss > $rootname.horiz

if ($status != 0) then
    echo "FATAL: Error whilst running psipass2 - script terminated!"
    exit $status
endif

# Remove temporary files

echo Cleaning up ...
\rm -f $tmproot.* error.log

echo "Final output files:" $rootname.ss2 $rootname.horiz
echo "Finished."