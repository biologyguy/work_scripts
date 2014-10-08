#!/usr/bin/python

from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
import os
from subprocess import Popen, PIPE, STDOUT
import re
import argparse
import get_common_ancestor

parser = argparse.ArgumentParser(prog="ancRecon2real.py", description="")
parser.add_argument('-p', '--pagan_dir', help='Path to the folder where pagan outfiles are stored. Don\'t put anything else in here...')
parser.add_argument('-s', '--subgroup_dir', help='Path to the folder where matrix_parse.py saved the subgroup files')
parser.add_argument('-o', '--outfile', help='Path to the new output file', default='ancRecon2real_output.txt')
incoming_args = parser.parse_args()


def get_common_ancestor(infile, anc_dir):
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
        
        for seq_id in sequence_ids_list:
            description = current_gen_seqs[seq_id].description.split(" ")[2]
            
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
            return SeqRecord(Seq(str(final_gen[parent_list[0]].seq), IUPAC.protein), id="Gen" + str(parental_gen).zfill(3) + ":Id" + final_gen[parent_list[0]].id + ":" + infile)


#alignment is a AlignIO.read() object
def calcIdentity(alignment):  
    j = 0  # counts positions in first sequence
    i = 0  # counts identity hits
    
    first_seq = alignment[0].seq
    for amino_acid in alignment[1].seq:
        if amino_acid == '-':
            pass
        else:
            if amino_acid == first_seq[j]:
                i += 1
        j += 1
    
    gap_strip = str(first_seq).replace('-', '')
    percent = 100 * i / float(len(gap_strip))
    return percent


def AlignmentIdentity(seq1, seq2):  # seq1 and seq2 are SeqRecord objects
    mafft_temp = open("out1.tmp", "w")
    mafft_temp.write(str(">" + seq1.id + "\n" + seq1.seq + "\n>" + seq2.id + "\n" + seq2.seq))
    mafft_temp.close()
    
    event = Popen('einsi --quiet --thread 20 out1.tmp > out2.tmp', shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    event.communicate()
    alignment = AlignIO.read("out2.tmp", "fasta")
    
    os.remove("out1.tmp")
    os.remove("out2.tmp")
    return calcIdentity(alignment)

path, dirs, num_alignement_groups = os.walk(incoming_args.pagan_dir).next()  # get the list of files from the pagan directory

num_files = len(num_alignement_groups)
num_files = (num_files - 1) / 3  # The -1 is to get rid of ./, and div by 3 because 3 files are output by pagan --output-ancestors

file_counter = 1
final_output = open(incoming_args.outfile, "w", 1)
final_output.write("Subgroup#\t#Seqs_in_SG\t%ID_w/i/SG\t%ID_btw_anc&recon\tAve%ID_btw_anc&SG\tAve%ID_btw_recon&SG\n")

print(num_files)
sys.exit

while num_files > 0:
    # caluclate average %identity within the subgroup
    try:  # this is to allow for the file numbers to not be all sequential
        with open(incoming_args.pagan_dir + "subgroup" + str(file_counter) + ".fas"):
            pass
    except IOError:
        file_counter += 1
        continue
       
    print("Processing Subgroup " + str(file_counter))
    running_sum = 0.0
    subgroup = list(SeqIO.parse(incoming_args.subgroup_dir + "gen100_sub" + str(file_counter) + ".fasta", "fasta"))
    num_sequences = len(subgroup)
    num_alignments = 1.0
    while len(subgroup) > 0:
        next_record = subgroup.pop()
        for next_subgroup in subgroup:  # recursively calling the ever shrinking subgroup list
            sys.stdout.write("\rWithin Subgroup Alignment " + str(int(num_alignments)) + " of " + str((num_sequences - 1) * ((num_sequences - 1) + 1) / 2) + " (" + str(num_sequences) + " sequences)")
            sys.stdout.flush()
            running_sum += AlignmentIdentity(next_record, next_subgroup)
            num_alignments += 1
    print("")
    average_ident = running_sum / num_alignments
    final_output.write("gen100_sub" + str(file_counter) + "\t" + "%.0f" % num_sequences + "\t" + "%.2f" % average_ident + "\t")
    
    #get the root ancestor of the subgroup as reconstructed by pagan as well as the real ancestor and determine their % ident
    pagan_tree_file = open(incoming_args.pagan_dir + "subgroup" + str(file_counter) + ".anctree", "r")
    ancestor_nodes = re.findall(r"#[0-9]+#", pagan_tree_file.next())
    pagan_tree_file.close()
    
    pagan_root_node = ancestor_nodes[len(ancestor_nodes) - 1]
    pagan_fasta = SeqIO.to_dict(SeqIO.parse(incoming_args.pagan_dir + "subgroup" + str(file_counter) + ".fas", "fasta"))
    seq = str(pagan_fasta[pagan_root_node].seq).replace("-", "")
    pagan_root_seq = SeqRecord(Seq(seq, IUPAC.protein), id=pagan_root_node)
    
    real_ancestor = get_common_ancestor(incoming_args.subgroup_dir + "gen100_sub" + str(file_counter) + ".fasta", incoming_args.subgroup_dir)
    
    alignment_ident = AlignmentIdentity(pagan_root_seq, real_ancestor)
    final_output.write("%.2f" % alignment_ident + "\t")
    
    #Calculate the average % ident between the real and Recon ancestors vs the sequences in the subgroup
    subgroup = list(SeqIO.parse(incoming_args.subgroup_dir + "gen100_sub" + str(file_counter) + ".fasta", "fasta"))  # need to recall this, because I destroyed the list earlier
    real_running_sum = 0.0
    recon_running_sum = 0.0
    counter = 1
    for next_subgroup in subgroup:
        sys.stdout.write("\rReal & Recon to Subgroup Alignment " + str(counter) + " of " + str(num_sequences))
        sys.stdout.flush()
        real_running_sum += AlignmentIdentity(next_subgroup, real_ancestor)
        recon_running_sum += AlignmentIdentity(next_subgroup, pagan_root_seq)
        counter += 1
    print("")
    real_average_ident = real_running_sum / num_sequences
    recon_average_ident = recon_running_sum / num_sequences
    
    final_output.write("%.2f" % real_average_ident + "\t%.2f" % recon_average_ident + "\n")
    file_counter += 1
    num_files -= 1
    

final_output.close()
print("\nDone. You can find your output at " + incoming_args.outfile)


#output = event.stdout.read()
#let x=1; for y in ./gen100_sub*; do echo gen100_sub$x.fasta\n; pagan --seqfile gen100_sub$x.fasta --outfile pagan/output$x --output-ancestors; let x++; done

