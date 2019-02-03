#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import division
from sys import argv


"""
Author:Xinyuan Min [950829573070]

Description: this is a script to implement global pairwise sequence alignment algorithm
"""
# import statements here


# functions between here and __main__

blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""


def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int, parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

# write your own functions below here
class Cell(object):
    """This module is to cretat Cell() class to store the relevant information related to a cell,
    such as the corresponding res1 & res2， similarity scores， the alignment path（cell origin）
    """
    def __init__(self, res1, res2, score, i = 0, j = 0, res1_res2 = 0, prepath = []):
        super(Cell, self).__init__()
        self.res1 = res1
        self.res2 = res2
        self.score = score
        self.res1_res2 = res1_res2
        self.prepath = prepath
        self.i = i
        self.j = j

    def calculateScore(self, seq1, seq2, diagonal, left, above, gap_penalty, end_penalty):
        """Return a integer which represents the optimal alignment score for each cell

        diagonal: Class object, the upper left cell
        left: Class object, the adjacent left cell
        above: Class object, the adjacent above cell
        gap_penalty: a positive integer which represents gap penalty
        end penalty： a positive integer which represents the end gap penalty

        This function is to calculate the optimal align score of a cell
        """
        score_list = []
        # consider end gap penalty
        if self.i == len(seq1) and self.j == len(seq2):
            score_list.append({'score': diagonal.score + self.res1_res2, 'origin': 'diagonal'})
            score_list.append({'score': left.score - end_penalty, 'origin': 'left'})
            score_list.append({'score': above.score - end_penalty, 'origin': 'above'})
        # when a gap is not at the end of alignment, only gap_penalty is used
        else:
            score_list.append({'score': diagonal.score + self.res1_res2, 'origin': 'diagonal'})
            score_list.append({'score': left.score - gap_penalty, 'origin': 'left'})
            score_list.append({'score': above.score - gap_penalty, 'origin': 'above'})            
        # sort the score_list in descending order
        score_list = sorted(score_list, key= lambda x: x['score'], reverse = True)
        # choose the highest score from the three options
        self.score = score_list[0]['score']
        # record the alignment origin of each cell
        self.prepath = []
        self.prepath.append(score_list[0]['origin'])
        # record multiple origins
        if score_list[1]['score'] == score_list[0]['score']:
            self.prepath.append(score_list[1]['origin'])
        if score_list[2]['score'] == score_list[0]['score']:
            self.prepath.append(score_list[2]['origin'])

    def getscore(self, res1, res2):
        """Return similarity score from BLOSUM62 matrix for two residues
    `
        res1: string, amino acid
        res2: string, amino acid

        This is a copy of score() function from the skeleton the teacher offered.
        I moved the socre() function under Class Cell() to assign the similarity score to a Class object 
        """
        lookup1 = BLOSUM62_ORDER[res1]
        lookup2 = BLOSUM62_ORDER[res2]
        self.res1_res2 = BLOSUM62_MATRIX[lookup1][lookup2]


def initiateMatrix(seq1, seq2, gap_penalty, end_penalty):
    """Return a list of list (initial matrix), the type of element in the inner list is Class object

    seq1: string, amino acid sequence
    seq2: string, amino acid sequence
    gap_penalty：a positive integer which represents gap penalty
    end penalty： a positive integer which represents the end gap penalty

    This function is to initiate an intial matrix 
    """    
    seq1 = '_' + seq1
    seq2 = '_' + seq2
    initial_matrix = []
    for idx1, res1 in enumerate(seq1):
        row = []
        initial_matrix.append(row)
        for idx2, res2 in enumerate(seq2):
            cell = Cell(res1, res2, 0, idx1, idx2)
            row.append(cell)
            # assign score attribute to each cell in the first row
            if idx1 == 0:
                cell.score = - end_penalty * idx2
            # assign score attribute to each cell in the first column
            if idx2 == 0:
                cell.score = - end_penalty * idx1
    return initial_matrix


def fillMatrix(seq1, seq2, gap_penalty, end_penalty):
    """Return a list of list (filled matrix), this matrix store the optimal score and alignment
    origin of each cell

    seq1: string, amino acid sequence
    seq2: string, amino acid sequence
    gap_penalty：a positive integer which represents gap penalty
    end penalty： a positive integer which represents the end gap penalty

    This function is to assign every cell (Class object) in the filled_matrix a optimal score,
    as well as to record the alignment origin of each cell
    """    
    filled_matrix = initiateMatrix(seq1, seq2, gap_penalty, end_penalty)
    seq1 = '_' + seq1
    seq2 = '_' + seq2
    for idx1, res1 in enumerate(seq1):
        for idx2, res2 in enumerate(seq2):
            cell = filled_matrix[idx1][idx2]
            cell.prepath = []
            if idx1 > 0 and idx2 > 0:
                # obtain the similarity score form blosum62 matrix
                cell.getscore(res1, res2)
                # get the optimal score among the three paths 
                cell.calculateScore(seq1, seq2, filled_matrix[idx1-1][idx2-1], filled_matrix[idx1][idx2-1],
                                    filled_matrix[idx1-1][idx2], gap_penalty, end_penalty)
    return filled_matrix


def printMatrix(matrix):
    """This function is to print a matrix.
    The type of element in the filled_matrix is CLass object, thus cannot be printed directly.
    By using this function, we can print the optimal socre in the filled_matrix easily
    """
    printed_matrix = [[[] for j in range(len(matrix[0]))] for i in range(len(matrix))]
    for idx_row, row in enumerate(matrix):
        for idx_col, cell in enumerate(row):
            printed_matrix[idx_row][idx_col] = cell.score
    return printed_matrix



def traceback(matrix):
    """Return two strings which represent the aligned seq1 and seq2

    matrix: list of list which represents the filled_matrix

    This function is to trace the direction of arrows and return the aligned sequences
    """
    aligned_seq1 = ''
    aligned_seq2 = ''
    # get the index of the bottom right cell in the filled_matrix
    i = len(matrix) - 1
    j = len(matrix[0]) - 1
    # len_1 & len_2 are used to record the number of residues that have been added 
    len_1 = 0
    len_2 = 0
    while i > 0 and j > 0:
        cell = matrix[i][j]
        # res1 & res2 matchs, add both residues to algined_seq respectively
        if cell.prepath == ['diagonal']:
            aligned_seq1 = cell.res1 + aligned_seq1
            aligned_seq2 = cell.res2 + aligned_seq2
            i -= 1
            j -=1
            len_1 += 1
            len_2 += 1
        # res2 aligned to a gap, add res2 to algined_seq2, add a gap to aligned_seq1
        elif cell.prepath == ['left']:
            aligned_seq1 = '_' + aligned_seq1
            aligned_seq2 = cell.res2 + aligned_seq2
            j -= 1
            len_2 += 1
        # res1 aligned to a gap, add res1 to algined_seq2, add a gap to aligned_seq2
        ## actually here ignored the situation of multiple origins, which is not correct
        else:
            aligned_seq1 = cell.res1 + aligned_seq1
            aligned_seq2 = '_' + aligned_seq2
            i -= 1
            len_1 += 1
    # add redundant residues or gaps to corresponding aligned sequences
    lack_1 = len(seq1) - len_1
    lack_2 = len(seq2) - len_2
    lack = lack_1 - lack_2
    aligned_seq1 = seq1[:max(lack, 0)] + '_' * max(-lack, 0) + aligned_seq1
    aligned_seq2 = seq2[:max(-lack, 0)] + '_' * max(lack, 0) + aligned_seq2

    return aligned_seq1, aligned_seq2

# Question 5
def percentIdentity(seq1, seq2, gap_penalty, end_penalty):
    """Return a 2-decimal percentage which represents the percentage identity of the two sequences

    seq1: string, amino acid sequence
    seq2: string, amino acid sequence
    gap_penalty：a positive integer which represents gap penalty
    end penalty： a positive integer which represents the end gap penalty

    This function is to get the percent identity of two aligned, this is
    calculted by multiplying the number of matches in the pair by 100, and
    dividing by the length of aligned sequence
    """
    matrix = fillMatrix(seq1, seq2, gap_penalty, end_penalty)
    aligned_seq1 = traceback(matrix)[0]
    aligned_seq2 = traceback(matrix)[1]
    identity = 0
    for idx in range(len(aligned_seq1)):
        # if the characters of the same position are identic, and are not gaps
        if aligned_seq1[idx] == aligned_seq2[idx] and aligned_seq1[idx] != '_':
            identity += 1
    decimal_identity = identity/len(aligned_seq1)
    percent_identity = "{:.2%}".format(decimal_identity)
    print(percent_identity)



if __name__ == "__main__":
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    # seq3: GPA1_ARATH
    seq3 = "MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA"
    # seq4: GPA1_ORYSI
    seq4 = "MGSSCSRSHSLSEAETTKNAKSADIDRRILQETKAEQHIHKLLLLGAGESGKSTIFKQIKLLFQTGFDEAELRSYTSVIHANVYQTIKILYEGAKELSQVESDSSKYVISPDNQEIGEKLSDIDGRLDYPLLNKELVLDVKRLWQDPAIQETYLRGSILQLPDCAQYFMENLDRLAEAGYVPTKEDVLYARVRTNGVVQIQFSPVGENKRGGEVYRLYDVGGQRNERRKWIHLFEGVNAVIFCAAISEYDQMLFEDETKNRMMETKELFDWVLKQRCFEKTSFILFLNKFDIFEKKIQKVPLSVCEWFKDYQPIAPGKQEVEHAYEFVKKKFEELYFQSSKPDRVDRVFKIYRTTALDQKLVKKTFKLIDESMRRSREGT"
    gap_penalty = int(argv[1])
    end_penalty = int(argv[2])                     
    # Question 3
    filled_matrix1 = fillMatrix(seq1, seq2, gap_penalty, end_penalty)
    print(printMatrix(filled_matrix1))
    print(traceback(filled_matrix1))
    # Question 4
    filled_matrix2 = fillMatrix(seq3, seq4, gap_penalty, end_penalty)
    print('alignment score of seq3 & seq4: ', printMatrix(filled_matrix2)[-1][-1])
    # Question 5
    percentIdentity(seq3, seq4, gap_penalty, end_penalty)



