#!/usr/bin/python3
# -*- coding: utf-8 -*- 
import subprocess
import os
import re
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

    
#set up arguments for the program
parser = argparse.ArgumentParser(prog="protein_family_search", description="Uses blastp to query whole predicted protein databases for a gene family. Requires a directory with'prot' blast database(s), a directory with the original fasta files that the blast databases were built from, and a fasta formated collection of curated example proteins from which to search (the more the better). Hit counts are stored and then output, so the higher the hit count, the more likely that the gene is of actual interest.")

parser.add_argument('reference_sequences', help='Provide a fasta file containing a collection of homologous sequences you want to use as a reference to pull out more sequences from the complete protein databases')
parser.add_argument('blastp_database', help='Path to directory with blastp databases. The databases need to be named as genus_species.')
parser.add_argument('original_fasta_file', help='path to directory with fasta files used to create the blastp databases. The genus and species names must be at the front of the file name, separated by an underscore (_).')
parser.add_argument('-n', '--num_threads', default='1', help='Slect the number of processing cores to use for the blastp search. Default is 1.')
parser.add_argument('-m', '--matrix', default='BLOSUM62', help='Specify the amino acid substitution matrix to use. Default is BLOSUM62.')
parser.add_argument('-e', '--evalue', default='0.01', help='Specify the e-value cutoff threashold. Default is 0.01')
parser.add_argument('-f', '--gene_family', help='Optional. Provide a gene family name that will be appended to the id of all hits.')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Where would you like the output files to go? Default is the working directory.', action='store')


incoming_args = parser.parse_args()

#initialize files
#incoming_args.output_dir = re.sub('/', '',incoming_args.output_dir)
strong_hits = open(incoming_args.output_dir + "/strong_hits.fasta", "w")
weak_hits = open(incoming_args.output_dir + "/weak_hits.fasta", "w")
strong_hits.close()
weak_hits.close()


#create a dictionary of SeqIO records using the accession # as key
def get_accession(record):
    parts = record.id.split(" ")
    return parts[0]


def perform_blast(blast_params):  # blast_params is a dictionary [orig_prot_db_fasta, blastp_db, query_file, evalue, num_threads, matrix, species]
    #Read query database into a dictionary so the records are available for output
    database_file = open(blast_params['orig_prot_db_fasta'])
    db_seq_dict = SeqIO.to_dict(SeqIO.parse(database_file, "fasta"), key_function=get_accession)
    database_file.close()

    blast_cline = NcbiblastpCommandline(db=blast_params['blastp_db'], query=blast_params['query_file'], evalue=blast_params['evalue'], outfmt=5, num_threads=blast_params['num_threads'], matrix=blast_params['matrix'])

    blast_result = blast_cline()

    xml_file = open(incoming_args.output_dir + "/temp_xml_file", "w+")
    xml_file.write(blast_result[0])
    xml_file.seek(0, 0)

    blast_iterator = NCBIXML.parse(xml_file)

    unique_idents_array = []

    #output_array_headings = ["gene_id", "#_hits", "ave_%id", "ave_E-value", "ave_bit_score"]
    output_array = []

    for record in blast_iterator:
        alignments = record.alignments[:]
    
        for next_alignment in alignments:
            title_info = next_alignment.title.split(' ')
            hsps = next_alignment.hsps[0]
        
            #skip result if it's a small match
            if hsps.query < 100:
                continue
            
            if title_info[1] in unique_idents_array:
                index = unique_idents_array.index(title_info[1])
                output_array[index][1] += 1
                output_array[index][2] += float(hsps.identities / len(hsps.query))
                output_array[index][3] += float(hsps.expect)
                output_array[index][4] += float(hsps.bits)
            
            else:            
                append_element = [title_info[1], 1, float(hsps.identities / len(hsps.query)), float(hsps.expect), float(hsps.bits)]
                output_array.append(append_element)
                unique_idents_array.append(title_info[1])

    final_output = blast_params['species'] + "\ngene_id\t#_hits\tave_%id\tave_E-value\tave_bit_score\n"

    strong_hits = open(incoming_args.output_dir + "/strong_hits.fasta", "a")
    weak_hits = open(incoming_args.output_dir + "/weak_hits.fasta", "a")

    for next_line in output_array:
        #set a threshold cutoff, so we don't get lots of low hit results    
        ave_id = str(round(next_line[2] / next_line[1], 3))
        ave_e_value = str(round(next_line[3] / next_line[1], 3))
        ave_bit_score = str(round(next_line[4] / next_line[1], 3))
        db_seq_dict[next_line[0]].id = blast_params['species'] + "_" + incoming_args.gene_family + "_" + db_seq_dict[next_line[0]].id
        
        if next_line[1] > 5:    
            final_output += next_line[0] + "\t" + str(next_line[1]) + "\t" + ave_id + "\t" + ave_e_value + "\t" + ave_bit_score + "\n"
            SeqIO.write(db_seq_dict[next_line[0]], strong_hits, 'fasta')
            
        else:
            SeqIO.write(db_seq_dict[next_line[0]], weak_hits, 'fasta')   

    strong_hits.close()
    weak_hits.close()

    final_output += "\n"   
    print(final_output)
    summary_file = open(incoming_args.output_dir + "/summary_output.txt", "w")
    summary_file.write(final_output)
    summary_file.close()
    
    #house keeping... get rid of the temp file
    xml_file.close()
    os.remove(incoming_args.output_dir + "/temp_xml_file")

directory_list = os.listdir(incoming_args.original_fasta_file)

for i in directory_list:
    if i[0] != '.':
        species = i.split(".")
        
        params_dic = dict(orig_prot_db_fasta=incoming_args.original_fasta_file + i,
                          blastp_db=incoming_args.blastp_database + species[0],
                          query_file=incoming_args.reference_sequences, evalue=incoming_args.evalue,
                          num_threads=incoming_args.num_threads, matrix=incoming_args.matrix, species=species[0])
        
        #print params_dic
        perform_blast(params_dic)


'''    
#[orig_prot_db_fasta, blastp_db, query_file, evalue, num_threads, matrix, species]

print "Length of query: " + str(alignments[0].length)
print "Length of alignment: " + str(len(alignments[0].hsps[0].query))       
        
print "Score: " + str(alignments[0].hsps[0].score)
print "Bits: " + str(alignments[0].hsps[0].bits)
print "Expect: " + str(alignments[0].hsps[0].expect)
print "Num_alignments: " + str(alignments[0].hsps[0].num_alignments)
print "Identities: " + str(alignments[0].hsps[0].identities)
print alignments[0].hsps[0].positives
print alignments[0].hsps[0].gaps
print alignments[0].hsps[0].strand
print alignments[0].hsps[0].frame
print alignments[0].hsps[0].query
print alignments[0].hsps[0].query_start
print alignments[0].hsps[0].match
print alignments[0].hsps[0].sbjct
print alignments[0].hsps[0].sbjct_start
'''