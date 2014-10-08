#!/usr/bin/python3

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import argparse
import os
from random import randint

parser = argparse.ArgumentParser(prog="dup_seqs", description="")

parser.add_argument('input_file')
incoming_args = parser.parse_args()


sequences = SeqIO.parse(incoming_args.input_file, "fasta")
output = open("dups_removed.fasta", "w")

temp_dir = "temp_" + str(randint(100, 10000000000))
os.makedirs(temp_dir)

first = True
for i in sequences:
    print(i.id)
    
    if first:
        print("HELLO!")
        SeqIO.write(i, output, 'fasta')
        first = False
        
    seq_file = open(temp_dir + "/seq_file", "w")
    SeqIO.write(i, seq_file, 'fasta')
    seq_file.close()
    
    output.close()
    
    os.system("makeblastdb -in dups_removed.fasta -dbtype prot -out " + temp_dir + "/blast")
    
    output = open("dups_removed.fasta", "a")
    
    blast_cline = NcbiblastpCommandline(db=temp_dir + "/blast", query=temp_dir + "/seq_file", outfmt=5)
    
    blast_result = blast_cline()
    
    xml_file = open(temp_dir + "/xml_file", "w+")
    xml_file.write(blast_result[0])
    xml_file.seek(0, 0)

    blast_iterator = NCBIXML.parse(xml_file)
    
    record = blast_iterator.i()
    alignments = record.alignments[:]

    for alignment in alignments: 
        hsps = alignment.hsps[0]
        alignment_title = alignment.title.split(" ")
        if i.id == alignment_title[1]:
            continue
        
        percent_ident = float(hsps.identities) / float(alignment.length)
        
        if percent_ident > 0.95:
            print("Deleted!\n")
            break
        
        else:
            SeqIO.write(i, output, 'fasta')
            print("Kept!\n")
            break
    
    xml_file.close()

output.close()
os.system("rm -r " + temp_dir)