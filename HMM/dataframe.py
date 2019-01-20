#!/usr/bin/env python

"""
Author: 

Description: this is a script to ...
"""
# Import statements
from sys import argv
from random import random
import pandas as pd
from pandas import DataFrame

# Function definitions
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

# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,\
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,\
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,\
    'W':0.013, 'Y':0.034 } 


def sample(events):
    """Return a key from dict based on the probabilities 

    events: dict of {key: probability}, sum of probabilities should be 1.0. 
    """
    key_options = list(events.keys())
    cum = [0 for i in key_options]

    cum[0] = events[key_options[0]]
    for i in range(1,len(events)):
        cum[i] = cum[i-1] + events[key_options[i]]
    # Should not be necessary, but for safety
    cum[len(cum)-1] = 1.0

    ref_point = random()
    pick = None
    i = 0
    while not pick and (i < len(cum)):
        if ref_point < cum[i]:
            pick = key_options[i]
        i = i + 1
    return pick

def matches(df):
    columns = len(df.columns)
    rows = len(df)
    matches = []
    for i in range(columns):
        counts = df[i].value_counts().to_dict()
        gaps = counts.setdefault('-', 0) # count gaps in column
        if gaps/rows < 0.5: # if number of gaps less than half
            matches.append('*')
        else:
            matches.append('')
    #match_df = pd.DataFrame(matches,ignore_index = True)
    total_matches= matches.count('*')
    df.loc[rows] = matches
    return df, total_matches
   
                          
def count_unique_in_column(df):
    
    print('count():',df[0].count())
    print('nunique():',df[0].nunique())
    print('size:',df[0].size)
    columns = len(df.columns)
    rows = len(df)
    info = {}
    for i in range(columns):
        if df[i].loc[rows] == '*':
            counts = df[i].value_counts().to_dict() # count distinct value by column
            if '-' in counts: # delect counts of gap
                counts.pop('-')
        info[i] =counts
     
    print(info)
   
    
           
if __name__ == "__main__":

    # implement main code here
    infile = 'test.fasta'
    data_test = fastaParser(infile)
    data_large = fastaParser('test_large.fasta')
    
    
    df = DataFrame(data_test)
    print(df)
    count_unique_in_column(df)
    print(matches(df)[0])
    


    
