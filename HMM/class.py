#!/usr/bin/env python3

"""
Author: 

Description: this is a script to ...
"""
# Import statements
from sys import argv
from random import random
import pandas as pd
from pandas import DataFrame
from collections import Counter
import math
# Function definitions
# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,\
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,\
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,\
    'W':0.013, 'Y':0.034 }

def fastaParser(filename):
    """this is a fasta parser to read the input data and return it as a list of strings

    filename: file name of input data
    """
    seqs = []
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                seqs.append(list(line.rstrip()))
    return seqs

def get_column(df):
    """Return a list of list, each inner list are the symbols in one column"""
    column_list = []
    nr_columns = len(df.columns) # nr. of columns in dataframe
    for i in range(nr_columns):
        column = df[i].tolist()
        column_list.append(column)
    return column_list

def mark_column(column_list):
    """Return a list of list, each list is marked to indicate the state"""
    for column in column_list:
        # return a dict with symbol as key, occurrence as value
        occurrence = Counter(column)
        # if occurrence of gap is less than half
        if occurrence.setdefault('-', 0)/ len(column) < 0.5:
            column.append('*') # if in match state, mark * in end of column
        else:
            column.append('.') # if not match, mark . in end of column
    return column_list
            
class Column():
    def __init__(self, col_ind, all_columns, e_m = {}):
        self.col_ind = col_ind
        self.symbols = all_columns[self.col_ind-1]
        self.state = self.symbols[-1]
        self.e_m = e_m
        
    def emission_probs(self, all_columns):
        if self.state == '*': # if in match state
        # return a dict with symbol as key, occurrence as value
            occurrence = Counter(self.symbols[:-1])
            if '-' in occurrence: # delete count of gaps if '-' in occurrence_dictionary
                occurrence.pop('-')
            for amino_acid in pa.keys():
                if amino_acid not in occurrence:
                    # add pesudo count 1 for those didn't occur in the column
                    occurrence[amino_acid] = 1
                else:
                    occurrence[amino_acid] += 1

            denominator = sum(occurrence.values()) # total counts of all amino_acids        
            for amino_acid, count in occurrence.items():
                self.e_m[amino_acid] = count/denominator # probability = count/total counts
        else: # if insertion state, take random probability
            self.e_m = pa
        return self.e_m
    
def transition_probs(ind, all_columns):
    transitions = ['mm', 'mi', 'im', 'ii', 'md', 'dm', 'dd']
    trans_probs = {}
    for key in transitions:
        trans_probs[key] = 0

    col = Column(ind, all_columns) # current column 
    if col.state == '.': # if current column in insertion state
        trans_probs = None
    else:
        next_col = Column(ind+1, all_columns) # next column
        if next_col.state == '*':
            for i in range(len(col.symbols)-1): # iterate over all symbols in the column
                if col.symbols[i] != '-': # if symbol in match state 
                    if next_col.symbols[i] != '-': # next symbol also in match state
                        trans_probs['mm'] += 1
                    else: # next symbol in deletion state
                        trans_probs['md'] += 1
                else: # if symbol in deletion state
                    if next_col.symbols[i] != '-': # next symbol in match state
                        trans_probs['dm'] += 1
                    else: # next symbol also in deletion state 
                       trans_probs['dd'] += 1

        else: # next column in insertion state
            for i in range(len(col.symbols)-1):
                if next_col.symbols[i] != '-': # next symbol in insertion state
                    if col.symbols[i] == '-': # symbol in deletion state
                        trans_probs['di'] += 1
                    else: # symbol in match state
                        trans_probs['mi'] += 1
         
                    
                
                
                
        
    

        
                        
                    
        
    
        
  
    
                    
            
                
    

    return trans_probs

  

    

        
        
    

if __name__ == "__main__":

    # implement main code here
    infile = 'test.fasta'
    data_test = fastaParser(infile)
    file_large = 'test_large.fasta'
    data_large = fastaParser(file_large)

    df_test = DataFrame(data_test) # covert to dafaframe
    
    all_columns = get_column(df_test)
    marked_columns = mark_column(all_columns)
    #print(mark_column(all_columns))


    
    #column1 = Column(marked_columns[0])
    column1 = Column(1, all_columns)
    column4 = Column(4, all_columns)
##    print(column1.emission_probs(all_columns))
##    print(column4.emission_probs(all_columns))
 
    # Question 3
##    for ind in range(1,len(marked_columns)+1):
##        col_obj = Column(ind, all_columns)
##        print(ind, [ round(math.log(i,10),2) for i in col_obj.emission_probs(all_columns).values()])
##

    print('col1:',transition_probs(1, all_columns))
    print('col2:',transition_probs(2, all_columns))
    print('col3:',transition_probs(3, all_columns))
    print('col4:',transition_probs(4, all_columns))
    print('col5:',transition_probs(5, all_columns))
    print('col6:',transition_probs(6, all_columns))
    print('col7:',transition_probs(7, all_columns))
    print('col8:',transition_probs(8, all_columns))
    print('col9:',transition_probs(9, all_columns))

