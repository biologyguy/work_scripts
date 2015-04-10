#!/usr/bin/python
# -*- coding: utf-8 -*-
from Bio import SeqIO, SubsMat
from Bio.Align import MultipleSeqAlignment
# from Bio.SubsMat import MatrixInfo
from scipy import stats
import numpy
import math
import sys
import re
import collections

import imp
SubsMat = imp.load_source(
    'SubsMat', '/Users/bondsr/Documents/biopython/Bio/SubsMat/__init__.py')
MatrixInfo = imp.load_source(
    'MatrixInfo', '/Users/bondsr/Documents/biopython/Bio/SubsMat/MatrixInfo.py')


class PSSM:
    def __init__(self, alignment):
        self.name = "Unnamed"

        # default amino acid priors: Grabbed them from
        # http://www.tiem.utk.edu/~gross/bioed/webmodules/aminoacid.htm and
        # then altered the matrix to include gap and X -> 1/22 ≈ 0.0445, and
        # then change all other values by i-(i*0.045)
        self.priors_matrix = {"A": 0.067, "C": 0.03, "D": 0.054, "E": 0.053, "F": 0.036, "G": 0.067, "H": 0.026,
                              "I": 0.035, "K": 0.066, "L": 0.069,"M": 0.016, "N": 0.04, "P": 0.046, "Q": 0.034,
                              "R": 0.038, "S": 0.074, "T": 0.056, "V": 0.062, "W": 0.012, "Y": 0.03, "X": 0.0445,
                              "-": 0.0445}

        # priors_matrix = {"A":0.074,"C":0.033,"D":0.059,"E":0.058,"F":0.04,"G":0.074,"H":0.029,"I":0.038,"K":0.072,
        # "L":0.076,"M":0.018,"N":0.044,"P":0.05,"Q":0.037,"R":0.042,"S":0.081,"T":0.062,"V":0.068,"W":0.013,"Y":0.033,
        # "X":0.001,"-":0.001}

        # arbitrarily set constant used to calculate prior corrected frequency
        # matrix (prevents any awkward 'divide by 0')
        # a value of 1.0 gives equal weight to observed frequency and prior
        # frequency; I don't know if this is good or not
        self.scaling_factor_k = 1.0

        # arbitrarily set constant used to scale the relative weight of
        # observed data versus the relative weight of the substitution matrix
        # (usually BLOSUM62, but PHAT for TM regions)
        self.scaling_factor_R = 2

        if isinstance(alignment, MultipleSeqAlignment):
            self.alignment = alignment

            # default substitution matrix is BLOSUM62 across entire sequence
            blosum62 = self.prep_seqmat(MatrixInfo.blosum62)
            # [[start,end,matrix],[start,end,matrix],...]
            self.assignMatrixPositions(
                [[0, len(self.alignment[0]) - 1, blosum62]])

        else:
            sys.exit(
                "You must provide a Bio.Align.MultipleSeqAlignment object of the sequences you want included in the "
                "PSSM. Use AlignIO.read (i.e., 'from Bio import AlignIO')")

    def build_pssm(self):
        # build pseudo-phylogentic weighting for all included sequences
        self.seq_weight = self._sequence_weights(self.alignment)

        # calculate raw frequency matrix f(i,j)
        self.freq_matrix = self._observed_residue_freqs(
            self.alignment, self.seq_weight)

        # Correct frequency matrix with prior frequencies f'(i,j)
        self.cor_freq_matrix = self._prior_correction(
            self.freq_matrix, self.priors_matrix)

        # LoD transform and apply SubMatrix f''(i,j)
        self.lod_matrix = self._lod_substitution_transform()

        return

    def print_pssm(self):
        try:
            output = "\t1\t2\t3\t4\t5\t6\t7\t8\n"
            for _next in self.lod_matrix:
                output += "%s\t" % _next
                for i in range(8):
                    output += "%\t" % str(round(self.lod_matrix[_next][i], 2))
                output += "\n"
            return output
        except Exception:
            return ValueError.message

    # eg # [(0,90,blosum62),(91,103,phat),...] The matrices must be full 20x20
    # SubsMat objects
    def assignMatrixPositions(self, positions):
        # This function doesn't do that much except validate that the right
        # information has been put into the positions matrix
        if positions[-1][1] != len(self.alignment[0]) - 1:
            sys.exit("ERROR: your list of lists at assignMatrixPositions did not end with the same length as the "
                     "alignement.\n%s\n%s" % (positions[-1][1], self.alignment[0]))

        start = 0
        self.matrix_assignment = []  # empty the current value
        for _next in positions:
            # Do some checking to make sure the data is formated correctly
            if len(_next) != 3:
                sys.exit("ERROR: provide a list of lists for positions at assignMatrixPositions()--> "
                         "[[start,end,matrix],[start,end,matrix],...]")

            if _next[0] != start:
                sys.exit("ERROR: the list of lists for positions at assignMatrixPositions() do not continuously span "
                         "the alignment.\nnext[0]: %s\tstart: %s" % (_next[0], start))

            if not isinstance(_next[2], SubsMat.SeqMat):
                sys.exit(
                    "ERROR: you must provide SeqMat matrix objects as the third parameter at assignMatrixPositions()")

            start = _next[1] + 1

            self.matrix_assignment.append(_next)
        return

    @staticmethod
    def prep_seqmat(matrix_info_obj):
        # Make a SeqMat object out of a matrix dictionary (like those in
        # MatrixInfo, strip out non-standard residues (e.g., B and Z) and add
        # in gap/missing penalties.
        matrix = SubsMat.SeqMat(matrix_info_obj, non_synon=True)
        _temp = SubsMat.SeqMat(matrix_info_obj, non_synon=True)

        kill_list = ["B", "J", "O", "U", "Z"]
        for _next in kill_list:
            matrix.alphabet.letters = matrix.alphabet.letters.replace(_next, "")
            if _next in matrix.ab_list:
                matrix.ab_list.remove(_next)

        matrix.alphabet.letters += "X-"
        for _next in matrix.alphabet.letters:
            matrix[(_next, "-")] = -4
            matrix[("-", _next)] = -4
            matrix[(_next, "X")] = -1
            matrix[("X", _next)] = -1

        matrix[("-", "-")] = 0
        matrix[("X", "X")] = 0

        for key in _temp:
            if key[0] in kill_list or key[1] in kill_list:
                del matrix[key]

        return matrix

    # topology file needs to be in fasta format, with all sequences
    def alignment_membranes(self, top_file, solu_matrix=MatrixInfo.blosum62, mem_matrix=MatrixInfo.slim161):
        # scans the alignment file, and compares it to topology. Pulling out the consensus membrane positions and
        # pushing them into a list. It also considers cases where a 'domain' is 1 or 2 residues, and in
        # those cases fuses the domains flanking the short one into a single domain.
        with open(top_file) as file:
            tops = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        tops_dic = {}
        for _next in tops:
            tops_dic[_next] = [list(tops[_next].seq), 0]

        matrix_positions = {"out": [], "mem": []}
        half_of_seqs = len(self.alignment) / 2
        for next_pos in range(len(self.alignment[0].seq)):
            memb_count = 0
            for next_align in self.alignment:
                if next_align.seq[next_pos] != "-":
                    if tops_dic[next_align.id][0][tops_dic[next_align.id][1]] == "M":
                        memb_count += 1
                    tops_dic[next_align.id][1] += 1

            if memb_count >= half_of_seqs:
                if len(matrix_positions['mem']) > 1:
                    if matrix_positions['mem'][-1] == next_pos - 2:
                        matrix_positions['mem'].append(
                            matrix_positions['out'].pop(-1))
                        matrix_positions['mem'].append(next_pos)

                    elif matrix_positions['mem'][-1] == next_pos - 3:
                        matrix_positions['mem'].append(
                            matrix_positions['out'].pop(-1))
                        matrix_positions['mem'].append(
                            matrix_positions['out'].pop(-1))
                        matrix_positions['mem'].append(next_pos)

                    else:
                        matrix_positions['mem'].append(next_pos)
                else:
                    matrix_positions['mem'].append(next_pos)

            else:
                if len(matrix_positions['out']) > 1:
                    if matrix_positions['out'][-1] == next_pos - 2:
                        matrix_positions['out'].append(
                            matrix_positions['mem'].pop(-1))
                        matrix_positions['out'].append(next_pos)

                    elif matrix_positions['out'][-1] == next_pos - 3:
                        matrix_positions['out'].append(
                            matrix_positions['mem'].pop(-1))
                        matrix_positions['out'].append(
                            matrix_positions['mem'].pop(-1))
                        matrix_positions['out'].append(next_pos)

                    else:
                        matrix_positions['out'].append(next_pos)
                else:
                    matrix_positions['out'].append(next_pos)

        final_array = []
        # first while loop builds the whole array except the very last domain
        # Note... If the space between 2 domains is <= 2 residues, merge into a
        # single domain
        # temp_matrix_positions = {"mem": [], "out": []}
        count = 1
        for _ in matrix_positions['mem']:
            if count > len(matrix_positions['mem']):
                break

            count += 1

        # prepare SubsMat objects for matrices
        solu_matrix = self.prep_seqmat(solu_matrix)
        mem_matrix = self.prep_seqmat(mem_matrix)

        while len(matrix_positions['out']) > 0 and len(matrix_positions['mem']) > 0:
            if matrix_positions['out'][0] < matrix_positions['mem'][0]:
                final_array.append(
                    [matrix_positions['out'].pop(0), 0, solu_matrix])
                while len(matrix_positions['out']) > 0:
                    if matrix_positions['out'][0] < matrix_positions['mem'][0]:
                        final_array[-1][1] = matrix_positions['out'].pop(0)
                    else:
                        break
            else:
                final_array.append(
                    [matrix_positions['mem'].pop(0), 0, mem_matrix])
                while len(matrix_positions['mem']) > 0:
                    if matrix_positions['mem'][0] < matrix_positions['out'][0]:
                        final_array[-1][1] = matrix_positions['mem'].pop(0)
                    else:
                        break

        # need to fill in the final domain with this next step. I'm pretty sure
        # there's a way to integrate this with the above while loop, but the
        # logic is escaping me...
        if len(matrix_positions['out']) == 0:
            final_array.append(
                [matrix_positions['mem'].pop(0), matrix_positions['mem'].pop(-1), mem_matrix])

        else:
            final_array.append(
                [matrix_positions['out'].pop(0), matrix_positions['out'].pop(-1), solu_matrix])

        self.assignMatrixPositions(final_array)

        return

    def _lod_substitution_transform(self):
        def apply_sub_matrix(position, aa):
            # apply the proper substitution matrix
            for _next in self.matrix_assignment:
                if _next[1] >= position:
                    sub_matrix = _next[2]
                    break

            output = 0.0
            for key in self.cor_freq_matrix:
                # Note:
                output += self.cor_freq_matrix[key][position] * \
                    sub_matrix[(aa, key)]
            return output

        lod_matrix = {"A": [], "C": [], "D": [], "E": [], "F": [], "G": [], "H": [], "I": [], "K": [], "L": [], "M": [],
                      "N": [], "P": [], "Q": [], "R": [], "S": [], "T": [], "V": [], "W": [], "Y": [], "X": [], "-": []}

        for next_aa in self.cor_freq_matrix:
            counter = 0
            for next_column in self.cor_freq_matrix[next_aa]:
                freq_by_prior = next_column / self.priors_matrix[next_aa]
                # freq_by_prior = next_column
                new_score = self.scaling_factor_R * \
                    math.log(freq_by_prior, 2) + \
                    (apply_sub_matrix(counter, next_aa))
                # new_score = [math.log(freq_by_prior,2),ApplySubMatrix(counter,next_aa)]
                lod_matrix[next_aa].append(round(new_score, 4))
                counter += 1

        return lod_matrix

    @staticmethod
    def _observed_residue_freqs(aln, seq_weight):
        freq_matrix = {"A": [], "C": [], "D": [], "E": [], "F": [], "G": [], "H": [], "I": [], "K": [], "L": [
        ], "M": [], "N": [], "P": [], "Q": [], "R": [], "S": [], "T": [], "V": [], "W": [], "Y": [], "X": [], "-": []}
        aln_len = len(aln[1].seq)

        for column in range(aln_len):
            aa_tally = {"A": 0.0, "C": 0.0, "D": 0.0, "E": 0.0, "F": 0.0, "G": 0.0, "H": 0.0, "I": 0.0, "K": 0.0,
                        "L": 0.0, "M": 0.0, "N": 0.0, "P": 0.0, "Q": 0.0, "R": 0.0, "S": 0.0, "T": 0.0, "V": 0.0,
                        "W": 0.0, "Y": 0.0, "X": 0.0, "-": 0.0}

            counter = 0
            for _next in aln:
                sequence = _next.seq
                if not re.match('[ACDEFGHIKLMNPQRSTUVWYX-]', sequence[column]):
                    aa_tally["X"] += seq_weight[counter]

                else:
                    aa_tally[sequence[column]] += seq_weight[counter]

                counter += 1

            for residue in aa_tally:
                freq_matrix[residue].append(aa_tally[residue])

        return freq_matrix

    def _prior_correction(self, freq_matrix, priors_matrix):
        corrected_matrix = {"A": [], "C": [], "D": [], "E": [], "F": [], "G": [], "H": [], "I": [], "K": [], "L": [
        ], "M": [], "N": [], "P": [], "Q": [], "R": [], "S": [], "T": [], "V": [], "W": [], "Y": [], "X": [], "-": []}

        for next_aa in freq_matrix:
            for next_column in freq_matrix[next_aa]:
                # The 1.0 in the denominator of this equation is the sum of
                # observed frequencies, which equals 1 instead of the number of
                # sequences in the alignment because of the phylogenetic
                # weighting
                new_freq = (
                    next_column + priors_matrix[next_aa] * self.scaling_factor_k) / (1.0 + self.scaling_factor_k)
                corrected_matrix[next_aa].append(new_freq)

        return corrected_matrix

    def write(self, path):
        seq_list = list(self.alignment[0].seq)
        with open(path, "w") as pssm_file:
            pssm_file.write("#PSSM for %s\n" % self.name)
            for _next in self.lod_matrix:
                if _next == "-" or _next == "X":
                    continue
                pssm_file.write(" " + _next)

            counter = 1
            for position in seq_list:
                pssm_file.write("\n" + str(counter) + " " + str(position))
                for _next in self.lod_matrix:
                    if _next == "-" or _next == "X":
                        continue
                    pssm_file.write(
                        " %.02f" % (self.lod_matrix[_next][counter - 1]))

                counter += 1
        return

    @staticmethod
    def _sequence_weights(aln):
        """ Original code written by Eric Talevich. http://permalink.gmane.org/gmane.comp.python.bio.devel/9522
        Weight aligned sequences to emphasize more divergent members.

        Returns a list of floating-point numbers between 0 and 1, corresponding to
        the proportional weight of each sequence in the alignment. The first list
        is the weight of the first sequence in the alignment, and so on. Weights
        sum to 1.0.

        Method: At each column position, award each different residue an equal
        share of the weight, and then divide that weight equally among the
        sequences sharing the same residue.  For each sequence, sum the
        contributions from each position to give a sequence weight.

        See Henikoff & Henikoff (1994): Position-based sequence weights.
        """
        def col_weight(column):
            """Represent the diversity at a position.

            Award each different residue an equal share of the weight, and then
            divide that weight equally among the sequences sharing the same
            residue.

            So, if in a position of a multiple alignment, r different residues
            are represented, a residue represented in only one sequence contributes
            a score of 1/r to that sequence, whereas a residue represented in s
            sequences contributes a score of 1/rs to each of the s sequences.
            """
            # Skip columns with all gaps or unique inserts
            if len([c for c in column if c not in '-.']) < 2:
                return [0] * len(column)
            # Count the number of occurrences of each residue type
            # (Treat gaps as a separate, 21st character)
            counts = collections.Counter(column)
            # Get residue weights: 1/rs, where
            # r = nb. residue types, s = count of a particular residue type
            n_residues = len(counts)    # r
            freqs = dict((aa, 1.0 / (n_residues * count))
                         for aa, count in iter(counts.items()))
            weights = [freqs[aa] for aa in column]
            return weights

        seq_weights = [0] * len(aln)
        col_weights = map(col_weight, zip(*aln))
        # Sum the contributions from each position along each sequence -> total
        # weight
        for col in col_weights:
            for idx, row_val in enumerate(col):
                seq_weights[idx] += row_val
        # Normalize
        scale = 1.0 / sum(seq_weights)
        seq_weights = [scale * wt for wt in seq_weights]
        return seq_weights


# This was interesting, but not as sophisticated as what I was later able
# to achieve with AlignMe.
class Compare:
    def __init__(self, pssm1, pssm2):
        self.max_scores_matrix = []
        self.window_size = 5
        self.position_penalty = 1

        self.pssm1 = pssm1
        self.pssm2 = pssm2

        self.pssm1_lod = pssm1.lod_matrix
        self.pssm2_lod = pssm2.lod_matrix

        # arbitrarily grab one of the dict elements for size
        self.pssm1_size = len(self.pssm1_lod['-'])
        self.pssm2_size = len(self.pssm2_lod['-'])
        self.linear_equation = "f(x)=x*" + str(float(self.pssm2_size) / float(self.pssm1_size)) + \
            "; set yrange [0:" + str(self.pssm2_size) + \
            "]; set xrange [0:" + str(self.pssm1_size) + "]"

    def compare(self, pssm1_positions=[1, -1], pssm2_positions=[1, -1]):
        def sort_scores(values, index):
            # sort the top score lists highest to lowest
            top_value_sort = []
            top_index_sort = []
            for x in range(len(values)):
                max_index = values.index(max(values))
                top_value_sort.append(values.pop(max_index))
                top_index_sort.append(index.pop(max_index))

            output = (top_value_sort, top_index_sort)
            return output

        for counter1 in range(self.pssm1_size - (self.window_size - 1)):
            top_score_values = []
            top_score_index = []

            for counter2 in range(self.pssm2_size - (self.window_size - 1)):
                # calculate unweighted score for pssm1 window j at pssm2 window
                # k
                current_score = 0.0
                for next_aa in self.pssm2_lod:
                    # strip out the gap and X scores, which are just inflating scores meaninglessly, and (i think) reducing the signal-to-noise
                    # if next_aa in ["-","X"]:
                    #    continue

                    for _next in range(self.window_size):
                        if self.pssm1_lod[next_aa][counter1 + _next] < 0 and self.pssm2_lod[next_aa][counter2 + _next] < 0:
                            continue
                        current_score += self.pssm1_lod[next_aa][
                            counter1 + _next] * self.pssm2_lod[next_aa][counter2 + _next]

                # weight score by the relative distance apart of j and k. This
                # assumes that the proteins are of somewhat equal length, and
                # homology exists across the majority of that length
                positional_weighting = abs(float(counter1) / float(self.pssm1_size) - float(
                    counter2) / float(self.pssm2_size)) * self.position_penalty
                current_score = current_score - \
                    abs(current_score * positional_weighting)

                top_score_values.append(round(current_score))
                top_score_index.append(counter2)

            sorted_scores = sort_scores(top_score_values, top_score_index)
            self.max_scores_matrix.append(sorted_scores)
        return

    def plot_data(self, num_cols):
        output = "Col#"

        for _next in range(num_cols):
            output += "\tscore" + str(_next + 1) + "\tindex" + str(_next + 1)

        output += "\n"

        counter = 1
        for _next in self.max_scores_matrix:
            output += str(counter)
            for x in range(num_cols):
                output += "\t" + str(_next[0][x]) + "\t" + str(_next[1][x])

            output += "\n"
            counter += 1
        return output

    def stats(self, column):
        score_list = []
        index_list = []
        line_index_list = []
        num_above = 0

        for _next in self.max_scores_matrix:
            score_list.append(_next[0][column - 1])
            index_list.append(_next[1][column - 1])
            line_index_list.append(
                float(self.pssm2_size) / float(self.pssm1_size) * float(_next[1][column - 1]))
            if _next[0][column - 1] > 600:
                num_above += 1

        average = numpy.mean(score_list)
        stdev = numpy.std(score_list)
        var = numpy.var(score_list)
        linear_regression = stats.linregress(index_list, score_list)

        output = {"average": average, "stdev": stdev, "var": var,
                  "linregress": linear_regression, "num_above": num_above}

        return output
