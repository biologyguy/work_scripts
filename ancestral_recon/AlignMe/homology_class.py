#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

# core python packages
import os
import time
import sys
import argparse
import re
import random
import types
import shutil
from multiprocessing import Process, Lock, Value, Array, cpu_count
from random import shuffle
from subprocess import Popen

# Numpy
import numpy as np

# BioPython
from Bio import SeqIO, AlignIO
# from Bio.SubsMat import MatrixInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

# my own custom class
from pssm import PSSM
from MyFuncs import run_multicore_function


class Homology():
    def __init__(self, fasta1, fasta2, out_dir):
        # instantiate a few variables
        self.trim_seqs1 = ""
        self.trim_seqs2 = ""
        self.top_seqs1 = ""
        self.top_seqs2 = ""
        self.einsi_seqs1 = ""
        self.einsi_seqs2 = ""
        self.alignme_transmembrane_weight = float()
        self.alignme_psipred_weight = float()
        self.alignme_pssm_weight = float()

        # get path info for input files
        self.fasta1 = os.path.abspath(fasta1)
        self.fasta2 = os.path.abspath(fasta2)        
                        
        # prep incoming sequence lists
        self.fasta1_dic = SeqIO.to_dict(SeqIO.parse(self.fasta1, 'fasta'))
        self.fasta1_dic = self._clean_fasta(self.fasta1_dic)
        self.fasta2_dic = SeqIO.to_dict(SeqIO.parse(self.fasta2, 'fasta'))
        self.fasta2_dic = self._clean_fasta(self.fasta2_dic)
        self.base_name1 = ".".join(str(self.fasta1.split("/")[-1]).split(".")[:-1])
        self.base_name2 = ".".join(str(self.fasta2.split("/")[-1]).split(".")[:-1])
        
        # Output directories  
        self.out_dir = os.path.abspath(out_dir)
        self.tmp_dir = "%s/temp" % self.out_dir
        self.pssm_dir = "%s/PSSM_FILES" % self.out_dir
        self.ss2_dir = "%s/SS2_FILES" % self.out_dir
        self.nnprf_dir = "%s/NN_PRF_FILES" % self.out_dir
        self.alignme_dir = "%s/ALIGNME_FILES" % self.out_dir
        self.shuffled_dir = "%s/SHUFFLED_FILES" % self.alignme_dir
        
        # Output files
        self.top_file1 = "%s/%s_top.fasta" % (self.out_dir, self.base_name1)
        self.top_file2 = "%s/%s_top.fasta" % (self.out_dir, self.base_name2)
        self.trim_seq_file1 = "%s/%s_trim.fasta" % (self.out_dir, self.base_name1)
        self.trim_seq_file2 = "%s/%s_trim.fasta" % (self.out_dir, self.base_name2)
        self.shuffle_scores_file = "%s/shuffle_scores.scores" % self.out_dir
        
        # multicore stuff
        self.stdout_lock = Lock()
        self.top_file_lock = Lock()
        self.trim_file_lock = Lock()
        self.shuffle_scores_file_lock = Lock()
        
        # initialize a few variables to 'None' for sanity
        self.pairwise_array = []
        self.shuffle_scores_dic = {}
        
    def _prepare_outfiles(self):
        # Create output directories
        self._mkdir(self.out_dir)
        self._mkdir(self.tmp_dir)
        self._mkdir(self.pssm_dir)
        self._mkdir(self.ss2_dir)
        self._mkdir(self.nnprf_dir)
        self._mkdir(self.alignme_dir)
        self._mkdir(self.shuffled_dir)
        
        # unlink top and trim files if they already exist...
        # The calling function opens file in append mode "a", so important to start with empty file
        self._clean_up([self.top_file1, self.top_file2, self.trim_seq_file1, self.trim_seq_file2])
    
    @staticmethod
    def _clean_fasta(fasta_dic):  # in case the incoming fasta file is an alignment
        for seq_id in fasta_dic:
            clean_seq = re.sub("-", "", str(fasta_dic[seq_id].seq))
            fasta_dic[seq_id].seq = Seq(clean_seq, generic_protein)
        
        return fasta_dic
    
    def execute_all(self):
        self._prepare_outfiles()
        
        # Run octopus
        print("Execute Octopus on %s:" % self.base_name1)
        run_multicore_function(self.fasta1_dic, self._octopus,
                               ["/usr/local/blastdbs/species_protein/Hydra_magnipapillata"])
        print("Execute Octopus on %s:" % self.base_name2)
        run_multicore_function(self.fasta2_dic, self._octopus,
                               ["/usr/local/blastdbs/species_protein/Hydra_magnipapillata"])
        self._clean_up(["%s/PSSM_PRF_FILES" % self.out_dir, "%s/RAW_PRF_FILES" % self.out_dir, 
                        "%s/exec_times-bloctopus.txt" % self.out_dir, "%s/exec_times-octopus.txt" % self.out_dir])
        
        # process top files
        print("%s top file processing:" % self.base_name1)
        run_multicore_function(self.fasta1_dic, self._process_top_files, [self.top_file1, self.trim_seq_file1])
        print("%s top file processing:" % self.base_name2)
        run_multicore_function(self.fasta2_dic, self._process_top_files, [self.top_file2, self.trim_seq_file2])
        
        # create MSAs with MAFFT
        print("MAFFT alignment of %s" % self.base_name1,
              Popen("einsi --thread -1 --quiet %s > %s/%s_einsi.fasta" % 
                    (self.trim_seq_file1, self.out_dir, self.base_name1), shell=True).wait())
        print(" --> DONE\n")
        print("MAFFT alignment of %s" % self.base_name2, Popen("einsi --thread -1 --quiet %s > %s/%s_einsi.fasta" % 
                                                               (self.trim_seq_file2, self.out_dir, self.base_name2), 
                                                               shell=True).wait())       
        print(" --> DONE\n")
        
        # Should I execute trimal at this point??  
        # prep a bunch of sequence dictionaries
        with open(self.trim_seq_file1, "r") as file:
            self.trim_seqs1 = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
            
        with open(self.trim_seq_file2, "r") as file:
            self.trim_seqs2 = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
        
        self.pairwise_array = []
        self.shuffle_scores_dic = {}
        for next1 in self.trim_seqs1:
            for next2 in self.trim_seqs2:
                self.pairwise_array.append((self.trim_seqs1[next1].id, self.trim_seqs2[next2].id))
                self.shuffle_scores_dic["%s-%s" % (self.trim_seqs1[next1].id, self.trim_seqs2[next2].id)] = False
        
        with open(self.top_file1, "r") as file:
            self.top_seqs1 = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
            
        with open(self.top_file2, "r") as file:
            self.top_seqs2 = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
         
        with open("%s/%s_einsi.fasta" % (self.out_dir, self.base_name1), "r") as seq_file:
            self.einsi_seqs1 = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
            
        with open("%s/%s_einsi.fasta" % (self.out_dir, self.base_name2), "r") as seq_file:
            self.einsi_seqs2 = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

        # generate PSSMs -- differentiates between transmembrane and non-TM regions for
        # BLOSUM62 and PHAT respectively (implemented in the PSSM class).
        align = AlignIO.read("%s/%s_einsi.fasta" % (self.out_dir, self.base_name1), "fasta")
        pssm1 = PSSM(align)
        pssm1.name = self.base_name1
        pssm1.alignmentMembranes(self.top_file1)
        pssm1.BuildPSSM()
        pssm1.write("%s/%s.pssm" % (self.out_dir, self.base_name1))
        
        align = AlignIO.read("%s/%s_einsi.fasta" % (self.out_dir, self.base_name2), "fasta")
        pssm2 = PSSM(align)
        pssm2.name = self.base_name2
        pssm2.alignmentMembranes(self.top_file2)
        pssm2.BuildPSSM()
        pssm2.write("%s/%s.pssm" % (self.out_dir, self.base_name2))
                
        # create custom pssms for each included sequence by deleting space rows
        # fasta1
        print("Custom PSSMs for %s" % self.base_name1)
        with open("%s/%s.pssm" % (self.out_dir, self.base_name1), "r") as pssm_file:
            pssm_lines = pssm_file.readlines()
            
        run_multicore_function(self.einsi_seqs1, self._pssm, [pssm_lines])
        
        # fasta2
        print("Custom PSSMs for %s" % self.base_name2)
        with open("%s/%s.pssm" % (self.out_dir, self.base_name2), "r") as pssm_file:
            pssm_lines = pssm_file.readlines()
            
        run_multicore_function(self.einsi_seqs2, self._pssm, [pssm_lines])

        # Now onto PSI-PRED
        print("Executing PSI-Pred on %s" % self.base_name1)
        # fasta1 first
        run_multicore_function(self.trim_seqs1, self._psi_pred)
        
        print("Executing PSI-Pred on %s" % self.base_name2)
        # fasta2 next
        run_multicore_function(self.trim_seqs2, self._psi_pred)

        # Run AlignMe on all pairwise combinations
        print("Pairwise AlignMe runs:")
        run_multicore_function(self.pairwise_array, self.alignme)
    
    def save_obj(self, save_file):
        if len(self.shuffle_scores_dic) == 0:
            print("Sorry, your can't save the shuffle_scores_dic until you've actually generated it!")
            return

        with open(os.path.abspath(save_file), "w") as file:
            file.write("#pairwise\tunshuffled\tshuffled\n")
            for next_score in self.shuffle_scores_dic:
                file.write("%s\t%s\t" % (next_score, self.shuffle_scores_dic[next_score][0]))
                np.savetxt(file, self.shuffle_scores_dic[next_score][1], newline='\t', fmt="%.2f")
                file.write("\n")

        print("%s saved." % save_file)

    def load_obj(self, save_file):
        if self.shuffle_scores_dic:
            print("Sorry, you've aleady populated shuffle_scores_dic in this object. You'll need to make a new one.")
            return
        try:       
            with open(os.path.abspath(save_file), "r") as file:
                self.shuffle_scores_dic = {}
                temp_ids = np.loadtxt(file, usecols=[0], dtype='str')
                file.seek(0, 0)
                temp_unshuf = np.loadtxt(file, usecols=[1])
                file.seek(0, 0)
                temp_shuf = np.loadtxt(file, usecols=np.arange(2, 102))
                
                counter = 0
                for next_id in temp_ids:
                    self.shuffle_scores_dic[next_id] = [round(temp_unshuf[counter], 1), np.array((temp_shuf[counter]))]
                    counter += 1
        except IOError:
            print("Oops, I can't find the .save file in your output directory. Have you previously saved "
                  "this particular project?")
          
    def get_alignment(self, pair):
        with open("%s/%s-%s.aln" % (self.alignme_dir, pair[0], pair[1]), "r") as file:
            alignment = AlignIO.read(file, 'clustal')
            return alignment 

    @staticmethod
    def _mkdir(path):
        if not os.path.exists(path):
            os.makedirs(path)
        return os.path.abspath(path)

    @staticmethod
    def _clean_up(dir_list):
        for next_dir in dir_list:
            if os.path.exists(next_dir):
                Popen("rm -R " + next_dir, shell=True).wait()
        return "Removed " + str(dir_list)
    
    def _octopus(self, seq_obj, blastdb):
        tmp_dir = self._mkdir(self.tmp_dir + "/" + seq_obj.id)
        
        with open("%s/%s.fa" % (tmp_dir, seq_obj.id), "w") as file:
            SeqIO.write(seq_obj, file, "fasta")
        
        with open(tmp_dir + "/NameFile.txt", "w") as file:
            seq_id = seq_obj.id.split(".")
            file.write(seq_id[0] + "\n")
        
        # with self.stdout_lock:
            print("Running bloctopus on %s, pid:%s" % (seq_obj.id, os.getpid()))
        Popen("bloctopus %s/NameFile.txt %s %s /usr/local/bin/blastall /usr/local/bin/blastpgp %s "
              "/usr/local/bin/makemat -P > /dev/null 2>&1" %
              (tmp_dir, tmp_dir, self.out_dir, blastdb[0]), shell=True).wait()
        
        # with self.stdout_lock:
        #    print("Running octopus on %s, pid:%s" % (seq_obj.id,os.getpid()) )            
        Popen("octopus %s/NameFile.txt %s/PSSM_PRF_FILES %s/RAW_PRF_FILES %s -N > /dev/null 2>&1" %
              (tmp_dir, self.out_dir, self.out_dir, self.out_dir), shell=True).wait()
        
        # do a little reformatting of the .nnprf files --> should probably re-implement this in python...
        # I just grabbed this perl line from the web
        Popen("perl -e 'local $/; $_ = <>; s/END(.*)//gs; print' %s/%s.nnprf > %s/%s.temp; mv %s/%s.temp %s/%s.nnprf" %
              (self.nnprf_dir, seq_obj.id, self.nnprf_dir, seq_obj.id, self.nnprf_dir, seq_obj.id, self.nnprf_dir,
               seq_obj.id), shell=True).wait()
            
        # self._clean_up([tmp_dir])
        return

    def _process_top_files(self, seq_obj, out_files):
        top_in_file = "%s/%s.top" % (self.out_dir, seq_obj.id)
        nnprf_file = "%s/%s.nnprf" % (self.nnprf_dir, seq_obj.id)
        
        with open(top_in_file, "r") as i_file:
            top_output = []
            file_conts = i_file.read()
            file_conts = file_conts.split('\n')
            file_conts.pop(0)
            file_conts = "".join(file_conts)
            membranes = re.finditer('[io]M+[io]', file_conts)

            top_output.append(file_conts)
            for membrane in membranes:
                tm = (membrane.start() + 1, membrane.end() - 1)
                top_output.append(tm)

        if len(top_output) != 5:
            # Kill this sequence if it doesn't have exactly four TM domains
            self._clean_up([nnprf_file])           
        
        else:
            # Only include the four transmembrane domains, and 20 residues after the final TM
            if len(top_output[0]) < top_output[4][1] + 20:
                end_point = len(top_output[0])
            else:
                end_point = top_output[4][1] + 20
            
            with self.top_file_lock:
                with open(out_files[0], "a") as top_file:
                    top_file.write(">%s\n%s\n\n" % (seq_obj.id, (top_output[0][0:end_point])))
            
            seq_obj.seq = seq_obj.seq[0:end_point]        
            
            with self.trim_file_lock:
                with open(out_files[1], "a") as seq_file:
                    SeqIO.write(seq_obj, seq_file, 'fasta')
                
            Popen("head -n %s %s > %s.tmp; mv %s.tmp %s;" %
                  (str(top_output[4][1] + 27), nnprf_file, nnprf_file, nnprf_file, nnprf_file), shell=True).wait()
        
        self._clean_up([top_in_file])
        return
    
    def _pssm(self, seq_obj, pssm_lines):
        with open("%s/%s.pssm" % (self.pssm_dir, seq_obj.id), "w") as pssm_out_file:
            pssm_out_file.write("\nTrimmed pssm file for %s\n%s" % (seq_obj.id, pssm_lines[0][1]))
            seq_list = list(seq_obj.seq)
            index = 2  # need to skip the first couple lines of the pssm file (header stuff)
            seq_position = 1
            for position in seq_list:
                if position != "-":
                    out_line = re.search("[A-Z\-].+", pssm_lines[0][index])
                    pssm_out_file.write("%s %s%s\n" % (str(seq_position), position, out_line.group(0)[1:]))
                    seq_position += 1                 
                index += 1
        return

    def _psi_pred(self, seq_obj):
        tmp_fasta = "%s/%s.fa" % (self.tmp_dir, seq_obj.id)
        with open(tmp_fasta, "w") as tmp_file:
            SeqIO.write(seq_obj, tmp_file, "fasta")
        
        Popen("psipred %s > /dev/null 2>&1" % tmp_fasta, shell=True).wait()
        self._clean_up([seq_obj.id + ".ss", seq_obj.id + ".horiz"])
        Popen("mv %s.ss2 %s/%s.ss2" % (seq_obj.id, self.ss2_dir, seq_obj.id), shell=True).wait()
        return
        
    def alignme(self, combination):
        sim_file_loc = "%s/%s-%s.simf" % (self.tmp_dir, combination[0], combination[1])
        
        self.alignme_transmembrane_weight = 4.2
        self.alignme_psipred_weight = 1.4
        self.alignme_pssm_weight = 0.2
        
        with open(sim_file_loc, "w") as sim_file:
            sim_file.write("weight: %s type: UniversalProfileSimilarity column: 3 headerlines: 7 profile1: %s/%s.nnprf "
                           "profile2: %s/%s.nnprf\n" % (self.alignme_transmembrane_weight,
                                                        self.nnprf_dir, combination[0], self.nnprf_dir, combination[1]))
            
            ss2_profile_text = " headerlines: 2 profile1: %s/%s.ss2 profile2: %s/%s.ss2\n" % \
                               (self.ss2_dir, combination[0], self.ss2_dir, combination[1])
            sim_file.write("weight: %s type: UniversalProfileSimilarity column: 4%s" %
                           (self.alignme_psipred_weight, ss2_profile_text))
            sim_file.write("weight: %s type: UniversalProfileSimilarity column: 5%s" %
                           (self.alignme_psipred_weight, ss2_profile_text))
            sim_file.write("weight: %s type: UniversalProfileSimilarity column: 6%s" %
                           (self.alignme_psipred_weight, ss2_profile_text))
            
            sim_file.write("weight: %s type: PositionSpecificSimilarity PSSM1: %s/%s.pssm PSSM2: %s/%s.pssm\n" %
                           (self.alignme_pssm_weight, self.pssm_dir, combination[0], self.pssm_dir, combination[1]))
        
        output_loc = "%s/%s-%s" % (self.alignme_dir, combination[0], combination[1])
        
        above_threshold_gap_opening_penalty = 5
        above_threshold_gap_extension_penalty = 3
        
        below_threshold_gap_opening_penalty = 3
        below_threshold_gap_extension_penalty = 1.5
        
        termini_gap_opening_penalty = 3
        termini_gap_extension_penalty = 1.5
        
        thresholds_for_penalties = 0.5
        
        strings = (self.tmp_dir, combination[0], self.tmp_dir, combination[1], sim_file_loc, output_loc, output_loc,
                   below_threshold_gap_opening_penalty, above_threshold_gap_opening_penalty,
                   below_threshold_gap_extension_penalty, above_threshold_gap_extension_penalty,
                   termini_gap_opening_penalty, termini_gap_extension_penalty, thresholds_for_penalties)
        
        Popen("alignme -fasta_file1 %s/%s.fa -fasta_file2 %s/%s.fa -similarity_score_file %s -output_aligned_profiles "
              "%s.prf -output_aligned_sequences %s.aln -below_threshold_gap_opening_penalty %s "
              "-above_threshold_gap_opening_penalty %s -below_threshold_gap_extension_penalty %s "
              "-above_threshold_gap_extension_penalty %s -termini_gap_opening_penalty %s "
              "-termini_gap_extension_penalty %s -thresholds_for_penalties %s" % strings, shell=True).wait()

        with open("%s.prf" % output_loc, "r") as handle:
            score = self.score_alignme(handle)
        self._clean_up(["%s.*" % output_loc])
        return score

    @staticmethod
    def score_alignme(alignme_file):
        file_lines = alignme_file.readlines()
        
        # clear out header rows
        while True:
            if not re.match("#", file_lines[0]):
                del file_lines[0]
            else:
                break   
        del file_lines[-1]
        
        tally = 0.0
        for _next in file_lines:
            regular = re.sub("\s+", ", ", _next)
            regular = re.sub("\?0", "0", regular)
            data = regular.split(", ")
            tally += abs(float(data[1]) - float(data[5])) + abs(float(data[2]) - float(data[6])) \
                + abs(float(data[3]) - float(data[7])) + abs(float(data[4]) - float(data[8]))
        
        return round(tally, 1)

    def make_shuffles(self, combination, num_shuffles=100):  # 'combination' is a tupel or list of names to be compared
        with open("%s/scores_file.dat" % self.shuffled_dir, "w") as scores_file:
            scores_file.write("# AlignMe scores for shuffled sequences\n#Seq1\tSeq2\tShuffle#\tScore")              
                   
        run_multicore_function(range(num_shuffles), self.shuffled_align, [combination[0], combination[1]])
        
        mv_dir = "%s/%s-%s" % (self.shuffled_dir, combination[0], combination[1])
        self._mkdir(mv_dir)
        Popen("mv %s/%s_*%s_*.prf %s" % (self.alignme_dir, combination[0], combination[1], mv_dir), shell=True).wait()
        Popen("mv %s/scores_file.dat %s" % (self.shuffled_dir, mv_dir), shell=True).wait()
        
        # in the shuffled_scores_dic, each index is a list => [real alignme score , nparray of all shuffled scores]
        with open("%s/%s-%s.prf" % (self.alignme_dir, combination[0], combination[1])) as orig_alignme:
            dict_id = "%s-%s" % (combination[0], combination[1])
            scores_file_data = np.loadtxt("%s/%s-%s/scores_file.dat" %
                                          (self.shuffled_dir, combination[0], combination[1]), usecols=[3])
            self.shuffle_scores_dic[dict_id] = [self.score_alignme(orig_alignme), scores_file_data]

    def make_all_shuffles(self, num_shuffles=100):
        # multi-core the actual shuffles, as oposed to the pairwise combinations
        counter = 1
        for _next in self.pairwise_array:
            print("  %s -- %s (%s of %s)" % (_next[0], _next[1], counter, len(self.pairwise_array)))
            self.make_shuffles(_next, num_shuffles)
            counter += 1
            # break
    
    def shuffled_align(self, iteration, args):  # args --> [id1,id2]
        # need to shuffle: - [sequence, pssm, nnprf, ss2]  Kill all temp files after process
        id1 = args[0]
        id2 = args[1]
        
        def get_top(topology_obj):
            membranes = re.finditer('[io]M+[io]', str(topology_obj.seq))
            output = [str(topology_obj.seq)]

            for _next in membranes:
                output += [(_next.start() + 1, _next.end() - 1)]
            return output
        
        top1 = get_top(self.top_seqs1[id1])
        top2 = get_top(self.top_seqs2[id2])
        
        def compile_profiles(seq_obj, _id):
            seq_pos_list = list(seq_obj.seq)
            output = list(range(len(seq_pos_list)))
            
            with open("%s/%s.pssm" % (self.pssm_dir, _id), "r") as pssm:
                pssm_list = pssm.readlines()
            
            with open("%s/%s.nnprf" % (self.nnprf_dir, _id), "r") as nnprf:
                nnprf_list = nnprf.readlines()
            
            with open("%s/%s.ss2" % (self.ss2_dir, _id), "r") as ss2:
                ss2_list = ss2.readlines()
            
            for _next in range(len(seq_pos_list)):
                output[_next] = [seq_pos_list[_next], nnprf_list[_next + 7], ss2_list[_next + 2], pssm_list[_next + 3]]
            
            return output
        
        seq1_pos_list = compile_profiles(self.trim_seqs1[id1], id1)
        seq2_pos_list = compile_profiles(self.trim_seqs2[id2], id2)
             
        def shuffle_seq_list(seq_pos_list, top):
            membranes = []
            non_membranes = []
            start = 0
            for _next in range(len(top) - 1):
                non_membranes += seq_pos_list[start:top[_next + 1][0]]
                membranes += seq_pos_list[top[_next + 1][0]:top[_next + 1][1]]
                start = top[_next + 1][1]
            non_membranes += seq_pos_list[start:]
               
            random.shuffle(membranes)
            random.shuffle(non_membranes)
            
            shuffled_seq = []
            
            non_start = 0
            non_end = 0
            mem_start = 0
            for _next in range(len(top) - 1):
                shuffled_seq += non_membranes[non_start:non_start + (top[_next + 1][0] - non_end)]
                shuffled_seq += membranes[mem_start:mem_start + (top[_next + 1][1] - top[_next + 1][0])]
                non_start += top[_next + 1][0] - non_end
                non_end = top[_next + 1][1]
                mem_start += top[_next + 1][1] - top[_next + 1][0]
            
            shuffled_seq += non_membranes[non_start:]
            
            return shuffled_seq
        
        shuffled1 = shuffle_seq_list(seq1_pos_list, top1)
        shuffled2 = shuffle_seq_list(seq2_pos_list, top2)
        
        shuf1_base_name = "%s_%04d" % (id1, iteration)
        shuf2_base_name = "%s_%04d" % (id2, iteration)
        
        def make_shuffled_files(_shuffle_seq_list, _id, base_name):
            with open("%s/%s.pssm" % (self.pssm_dir, base_name), "w") as pssm:
                pssm.write("\nTrimmed pssm file for shuffled%s\n  A C E D G F I H K M L N Q P S R T W V Y\n" % _id)
                counter = 0
                for _ in _shuffle_seq_list:
                    re_indexed = re.sub("[0-9]+", str(counter + 1), _shuffle_seq_list[counter][3], 1)
                    pssm.write(re_indexed)
                    counter += 1
            
            with open("%s/%s.nnprf" % (self.nnprf_dir, base_name), "w") as nnprf:
                nnprf.write("Sequence: %s\n\nLength of query sequence: %s\n\n\nSTART 1\nALPHA:   M       L       "
                            "G       I       -      <SPACE> <LABEL> <QUERY>\n" % (base_name, len(_shuffle_seq_list)))
                counter = 0
                for _ in _shuffle_seq_list:
                    re_indexed = re.sub("[0-9]+", str(counter + 1), _shuffle_seq_list[counter][1], 1)
                    nnprf.write(re_indexed)
                    counter += 1
            
            with open("%s/%s.ss2" % (self.ss2_dir, base_name), "w") as ss2:
                ss2.write("# PSIPRED via homology_class.shuffled_align.make_shuffled_files() \n\n")
                counter = 0
                for _ in _shuffle_seq_list:
                    re_indexed = '{:>31}'.format(re.sub("[0-9]+", str(counter + 1), _shuffle_seq_list[counter][2], 1).lstrip())
                    ss2.write(re_indexed)
                    counter += 1
            
            with open("%s/%s.fa" % (self.tmp_dir, base_name), "w") as fa:
                fa.write(">%s\n" % _id)
                seq = ""
                for _next in _shuffle_seq_list:
                    seq += _next[0]
                fa.write("%s\n" % seq)
            return
        
        # make tmp pssm
        make_shuffled_files(shuffled1, id1, shuf1_base_name)
        make_shuffled_files(shuffled2, id2, shuf2_base_name)

        score = self.alignme([shuf1_base_name, shuf2_base_name])
        
        with self.shuffle_scores_file_lock:
            with open("%s/scores_file.dat" % self.shuffled_dir, "a") as handle:
                handle.write("\n%s\t%s\t%s\t%s" % (id1, id2, iteration, score))
        
        self._clean_up(["%s/%s-%s" % (self.shuffled_dir, id1, id2), "%s/%s.pssm" % (self.pssm_dir, shuf1_base_name),
                        "%s/%s.pssm" % (self.pssm_dir, shuf2_base_name), "%s/%s.ss2" % (self.ss2_dir, shuf1_base_name),
                        "%s/%s.ss2" % (self.ss2_dir, shuf2_base_name), "%s/%s.nnprf" % (self.nnprf_dir, shuf1_base_name),
                        "%s/%s.nnprf" % (self.nnprf_dir, shuf2_base_name), "%s/%s.fa" % (self.tmp_dir, shuf1_base_name),
                        "%s/%s.fa" % (self.tmp_dir, shuf2_base_name)])
        return
