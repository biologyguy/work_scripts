#!/usr/bin/python3

import argparse
import sys
import re
from Bio import SeqIO


class TMHMMrecord:
    def __init__(self, tmhmm_output_line):
        
        tmhmm_data = tmhmm_output_line.split("\t")
        
        self.id = tmhmm_data[0]
        self.length = int(tmhmm_data[1][4:])
        self.pred_hel = int(tmhmm_data[4][8:])
        self.topology = []
        
        #topology data is in the form e.g., o48-70i124-146o227-249i294-316o
        #where the numbers represent the transmembrane region        
        topology_record = tmhmm_data[5][9:].split("-")
    
        #deal with sequences that contain tm domains
        if topology_record[0].strip() != 'o':
            first = topology_record.pop(0)        
            last = topology_record.pop().strip()
        
            #place the first record
            self.topology.append([first[0:1], 1, int(first[1:]) - 1])
        
            prev_record = first
            for next_span in topology_record:
                topo_nums = re.findall('[0-9]+', next_span)
                topo_letter = re.findall('[io]', next_span)
    
                self.topology.append(["t", int(prev_record[1:]), int(topo_nums[0])])
                self.topology.append([topo_letter[0], int(topo_nums[0]) + 1, int(topo_nums[1]) - 1])
            
                prev_record = topo_letter[0] + topo_nums[1]

            topo_nums = re.findall('[0-9]+', last)
            topo_letter = re.findall('[io]', last)
        
            self.topology.append(["t", int(prev_record[1:]), int(topo_nums[0])])
            self.topology.append([topo_letter[0], int(topo_nums[0]) + 1, self.length])
        
        #if no tm domains, set whole sequence 
        else:
            self.topology.append(["o", 1, self.length])

#Mus_musculus_Innexin_ENSMUSP00000124354    len=677    ExpAA=86.15    First60=9.46    PredHel=4    Topology=o48-70i124-146o227-249i294-316o    

#set up arguments for the program
#*************************************************#
parser = argparse.ArgumentParser(prog="truncate_transmembrane_proteins.py", description="Curates a collection of putative transmembrane proteins to ensure they meet specific requirements, and removes specified amounts of sequence from the ends of each protein.")

parser.add_argument('input_sequences', help='FASTA formated list of putative transmembrane proteins that was used to create TMHMM output file.')
parser.add_argument('tmhmm_input', help='This file should be the contents of the output created by uploading the input_sequences FASTA file to the TMHMM server (http://www.cbs.dtu.dk/services/TMHMM/) and selecting the \'One line per protein\' output format radio button. Do not copy the header info from the output page, only the actual results.')
parser.add_argument('output_file', help='Path to the file you want the results saved to.') 
parser.add_argument("-n", '--num_tm_domains', default=4, help='Choose the exact number of transmembrane domains that must be present in each sequence for it to be kept in the final output. Default = 4.', type=int)
parser.add_argument("-t", '--trim', default='CTERM', help='Specify what should be trimmed from the protein sequence. Trimming will occur from the terminal residue to the nearest residue of the nearest tm-domain. Which termini are trimmed can be selected with the options CTERM (default), NTERM, BOTH, or NONE.')

incoming_args = parser.parse_args()
#*************************************************#

possible_trim_types = {'CTERM', 'NTERM', 'BOTH', 'NONE'}

if incoming_args.trim.upper() not in possible_trim_types:
    sys.exit("Error. Unknown trim type '" + incoming_args.trim.upper() + "' supplied. Accepted options are CTERM, NTERM, BOTH, or NONE.")


input_sequences = SeqIO.parse(incoming_args.input_sequences, "fasta")
tmhmm_input = open(incoming_args.tmhmm_input, 'r')
output_file = open(incoming_args.output_file, 'w')

for next_seq in input_sequences:
    tmhmm_rec = TMHMMrecord(tmhmm_input.next())
    if tmhmm_rec.pred_hel != 4:
        continue
    
    if next_seq.id != tmhmm_rec.id:
        sys.exit("Error. There was a missmatch between the TMHMM file (ID: " + tmhmm_rec.id + ") and the FASTA file (ID: " + next_seq.id + "). Please ensure that all records are present in the exact same order in both files.")
        
    output = ">" + tmhmm_rec.id + "\n" + next_seq.seq[0:tmhmm_rec.topology[7][2]] + "\n"
    output_file.write(str(output))


tmhmm_input.close()
output_file.close()