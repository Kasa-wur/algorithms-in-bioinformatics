#!/usr/bin/env python

"""
Author: 

Description: this is a script to ...
"""
# Import statements
from sys import argv
from random import random

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
                seqs.append(line.rstrip())
    return seqs

def nr_match(data):
    total_seqs = len(data)  # count total nr. of sequences in data
    length = len(data[0])   # get the length of one sequence
    nr_match = 0
    for j in range(length):
        count_gap = 0
        for i in range(total_seqs):
            if data[i][j] == '-':
                count_gap += 1
        if count_gap/total_seqs < 0.5: # if nr. of gaps is less than half
            nr_match += 1               
    return nr_match

            
# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,\
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,\
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,\
    'W':0.013, 'Y':0.034 }            

class HMM():
    """HMM object to store an HMM model

    This object is designed to keep track of all HMM states, emissions, and 
    transitions. It may be used in your implementation, but may also be 
    ignored, and replaced by a data structure of choice
    """
    # Emission probabilities for the match and insert states
    e_m   = []; e_i   = pa; 
    
    # Transition probabilities from/to matches, inserts and deletions
    t_mm  = []; t_mi  = []; t_md = [];
    t_im  = []; t_ii  = []
    t_dm  = []; t_dd  = []; 
    
    def __init__(self,nmatches):
        """Initialize HMM object with number of match states
        
        nmatches: int, number of match states
        """
    
        self.nmatches = nmatches
        
        self.e_m   = [dict(pa) for i in range(0,nmatches)]
        for i in range(0,nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i   = pa;

        self.t_mm  = [0.0 for i in range(0,nmatches+1)]
        self.t_mi  = [0.0 for i in range(0,nmatches+1)]
        self.t_im  = [0.0 for i in range(0,nmatches+1)]
        self.t_ii  = [0.0 for i in range(0,nmatches+1)]
        self.t_md  = [0.0 for i in range(0,nmatches+1)]
        self.t_dm  = [0.0 for i in range(0,nmatches+1)]
        self.t_dd  = [0.0 for i in range(0,nmatches+1)]


        

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

    
           
if __name__ == "__main__":

    # implement main code here
    infile = 'test.fasta'
    data_test = fastaParser(infile)
    #for seq in data_test:
#        print(seq)

    data_large = fastaParser('test_large.fasta')
    for seq in data_large[:20]:
        print(seq)



    
