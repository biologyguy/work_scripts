#!/usr/bin/python
from Bio import SeqIO
import re
import sys


def GetCommonAncestors(infile,anc_dir):

    starting_sequences = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))
    sequence_ids_list = []

    for keys in starting_sequences:
        sequence_ids_list.append(keys)

    num_seqs = len(sequence_ids_list)

    incoming_gen = re.search("gen[0-9]{3}", infile)
    if incoming_gen is None:
        sys.exit("Could not determine the starting generation number from the input file provided")

    parental_gen = int(incoming_gen.group(0)[3:])

    while num_seqs > 1:
        parent_list = []
        
        if parental_gen == 0:
            first_gen = SeqIO.parse(anc_dir + "gen000.fasta", "fasta")
            sys.exit("Ancestral seq from Gen000: " + first_gen.next().seq)
        
        current_gen_seqs = SeqIO.to_dict(SeqIO.parse(anc_dir + "gen" + str(parental_gen).zfill(3) + ".fasta", "fasta"))
        
        for i in sequence_ids_list:
            description = current_gen_seqs[i].description.split(" ")[2]
            
            if description in parent_list:
                num_seqs -= 1
                continue
            
            parent_list.append(description)        
        
        sequence_ids_list = parent_list        
        parental_gen -= 1
                
        if len(parent_list) == 1:
            final_gen = SeqIO.to_dict(SeqIO.parse(anc_dir + "gen" + str(parental_gen).zfill(3) + ".fasta", "fasta"))
            return ">Gen" + str(parental_gen).zfill(3) + ":Id" + final_gen[parent_list[0]].id + ":" + infile + "\n" + final_gen[parent_list[0]].seq
