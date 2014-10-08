#!/usr/bin/python3
from Bio import SeqIO

sequences_handle = SeqIO.parse('strong_hits.fasta', 'fasta')
seq_dict = SeqIO.to_dict(sequences_handle)

remove_list = open('delete_claudins', 'r')

for next_line in remove_list:
    next_line = next_line.strip()
    del seq_dict[next_line]
    
output_file = open('strong_hits_trimmed.fasta', 'w')

for i in seq_dict:
    SeqIO.write(seq_dict[i], output_file, 'fasta')


remove_list.close()
output_file.close()