#!/usr/bin/python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
from subprocess import Popen, PIPE, STDOUT
import re
import os


def GetCommonAncestor(infile, anc_dir):
    sys.stdout.write("Retrieving common ancestor ")
    sys.stdout.flush()
    
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
            final_gen = list(SeqIO.parse(anc_dir + "gen000.fasta", "fasta"))
            print("")
            return SeqRecord(Seq(str(final_gen[0].seq), IUPAC.protein), id="Gen000:Id" + final_gen[0].id + ":" + infile)
                    
        current_gen_seqs = SeqIO.to_dict(SeqIO.parse(anc_dir + "gen" + str(parental_gen).zfill(3) + ".fasta", "fasta"))
        
        for i in sequence_ids_list:
            description = current_gen_seqs[i].description.split(" ")[2]
            
            if description in parent_list:
                num_seqs -= 1
                continue
            
            parent_list.append(description)        
        
        sequence_ids_list = parent_list        
        parental_gen -= 1
        
        number_parents = len(parent_list)
        sys.stdout.write("\rRetrieving common ancestor ---> " + str(number_parents).zfill(3) + " sequences left to converge, at Gen " + str(parental_gen).zfill(3))
        sys.stdout.flush()
        
        if number_parents == 1:
            final_gen = SeqIO.to_dict(SeqIO.parse(anc_dir + "gen" + str(parental_gen).zfill(3) + ".fasta", "fasta"))
            print("")
            return SeqRecord(Seq(str(final_gen[parent_list[0]].seq), IUPAC.protein), id="Gen" + str(parental_gen).zfill(3) + ":Id" + final_gen[parent_list[0]].id + ":" + infile.split("/")[-1])

recon_output_file = open("recon_sequences.fas", "w")
real_output_file = open("real_sequences.fas", "w")

path, dirs, num_alignement_groups = os.walk("seq_files/pagan/").next()  # get the list of files from the pagan directory

for next_file in num_alignement_groups:
    if next_file.find('fas') != -1:
        subgroup = SeqIO.to_dict(SeqIO.parse("seq_files/pagan/" + next_file, 'fasta'))
        recon_id = "#" + str((len(subgroup) - 1) / 2) + "#"
        root_recon = subgroup[recon_id]
        root_recon.id = re.search("subgroup[0-9]+", next_file).group(0) + ":" + root_recon.id
        clean_seq = ""
        for char in root_recon.seq:
            if char != "-":
                clean_seq += char
        clean_seq = Seq(clean_seq, IUPAC.protein)
        root_recon.seq = clean_seq    
        
        real_ancestor = GetCommonAncestor("seq_files/gen100_sub" + re.search("[0-9]+", next_file).group(0) + ".fasta", "seq_files/")
        
        print(root_recon, "\n", real_ancestor)
        
        SeqIO.write(root_recon, recon_output_file, "fasta")
        SeqIO.write(real_ancestor, real_output_file, "fasta")
        

recon_output_file.close()
real_output_file.close()
sys.exit("All done")
