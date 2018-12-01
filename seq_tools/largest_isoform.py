#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from buddysuite import SeqBuddy as Sb
from buddysuite import buddy_resources as br
import argparse


def main():

    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="largest_isoform", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mLargest Isoform\033[m
  Select only the largest isoform from Augustus protein models  

  Pass in a file containing protein sequences with the .t# suffix at the end
  of each sequence ID.
  
\033[1mUsage\033[m:
  largest_isoform.py "/path/to/sequences"
''')

    # Positional
    parser.add_argument("sequences", help="Specify a sequence file")
    parser.add_argument("-i", "--in_place", action="store_true",
                        help="Overwrite original file. Be sure you want to do this!!")
    in_args = parser.parse_args()

    final_records = []
    seqs = Sb.SeqBuddy(in_args.sequences)
    seqs = Sb.order_ids(seqs)

    iso_id = ".".join(seqs.records[0].id.split(".")[:-1])
    max_seq = seqs.records[0]
    for rec in seqs.records:
        cur_id = ".".join(rec.id.split(".")[:-1])
        if cur_id != iso_id:
            iso_id = cur_id
            final_records.append(max_seq)
            max_seq = rec
        else:
            if len(max_seq.seq) < len(rec.seq):
                max_seq = rec

    seqs.records = final_records
    if in_args.in_place:
        seqs.write(in_args.sequences)
    else:
        print(seqs)


if __name__ == '__main__':
    main()