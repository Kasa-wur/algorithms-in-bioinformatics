#!/usr/bin/env python3

"""
Author: Xinyuan Min 950829573070

Description: this is a script to build profile HMM model and generate random sequence
"""
# Import statements
from sys import argv
from random import random
from collections import Counter

# Function definitions
# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,\
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,\
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,\
    'W':0.013, 'Y':0.034 }

def fastaParser(filename):
    """Return list of strings， each string is a sequence
    
    filename： string， e.g.: "test.fasta"
    this is a fasta parser to read the input data and return it as a list of strings
    filename: file name of input data
    """
    seqs = []
    seq = ''
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                seq += '>'
            else:
                seq += line.rstrip('\n')
    seq = seq.lstrip('>')
    seqs = seq.split('>')
    return seqs


def mark_state(data):
    """Return list of list, [['V','V','V','F','I','*'],[...],..]

    data: list of strings, each string is a sequence

    each inner list are amino acids that share the same position, as well as
    the state of that position. '*':match state, '.' :insertion state
    """
    column_list = []
    for i in range(len(data[0])):
        column = []
        for row in data:
            column.append(row[i])
        # return a dict with symbol as key, occurrence as value
        occurrence = Counter(column)
        if occurrence.setdefault('-', 0)/ len(column) < 0.5:
            column.append('*') # if in match state, mark * in end of column
        else:
            column.append('.') # if not match, mark . in end of column
        column_list.append(column)
    return column_list

def emission_probs(all_columns):
    """Return dict of dict,{0: {'V':0.222,'F':0.074,...}, 1: {'V':..,'F':..},..}

    all_columns: list of list, each element is a list of string

    this function calculate the emission probabilities of each position/column
    stored in a dict of dict, each position is a key and emission probabilies is
    value, emission probabilies is again a dict, amino acid as key and corre-
    sponding probability as value.
    """
    all_em = {}
    for ind, column in enumerate(all_columns):
        emission_probs = {}
        if column[-1] == '*': # if column in match state
            # return a dict with symbol as key, occurrence as value
            occurrence = Counter(column[:-1])
            if '-' in occurrence: # delete count of gaps 
                occurrence.pop('-')
            for amino_acid in pa.keys():
                if amino_acid not in occurrence:
                    # add pesudo count 1 for those didn't occur in the column
                    occurrence[amino_acid] = 1
                else:
                    occurrence[amino_acid] += 1
            # total counts of all amino_acids                       
            denominator = sum(occurrence.values())
            for amino_acid, count in occurrence.items():
                emission_probs[amino_acid] = count/denominator 
        else: # if in insertion state, take random probability
            emission_probs = pa
        all_em[ind] = emission_probs
    return all_em

def transition_dict():
    """Return a dictionary, transition state as keys, and 0 as values

    initiate empty transition state count dict, i.e. {'mm':0, 'mi':0,
    'im':0,...,'di':0}
    """
    transitions = ['mm', 'mi', 'im', 'ii','id', 'md', 'dm', 'dd','di' ]
    trans_probs = {}
    for key in transitions:
        trans_probs[key] = 0
    return trans_probs
    
def get_states(all_columns, ind):
    """Return a list of string， indicate states of symbols of one column
    for example，['m','m','d','m'] means [match, match, deletion,match]
    """
    col = all_columns[ind]
    col_state = []
    if col[-1] == '*': # column in match state
        for symbol in col[:-1]: # for symbol in col
            if symbol != '-': 
                col_state.append('m')
            else:
                col_state.append('d')
    else: # column in insertion state
        for symbol in col[:-1]:
            if symbol != '-': 
                col_state.append('i') 
            else: # if there is insertion
                col_state.append('')             
    return col_state

        
def match_match(all_columns, ind):
    """Return a diction，transition（string） as key，and counts（int）as value

    all_columns:list of list, each element is a list of string
    ind: int, index of current postion/ column

    this is support function for calculating transition probabilities.
    this function is only appiled when current and next position/column
    are both in match state.
    """
    col = all_columns[ind]  # get current column
    next_col = all_columns[ind+1] # get next column
    col_state = get_states(all_columns, ind) # get state of current column
    next_col_state = get_states(all_columns, ind+1) # get next column state
    # combine two state together, e.g., 'm'+'d' --> 'md'
    transition = list(map(lambda x,y: x+y, col_state, next_col_state))
    counts = Counter(transition) # count occurrence of each transition
    trans_counts = transition_dict().copy() # copy empty transition dict
    trans_counts.update(counts) # update the counts of each transition
    return trans_counts

def match_insertion(all_columns, ind):
    """Return a diction，transition（string） as key，and counts（int）as value

    all_columns:list of list, each element is a list of string
    ind: int, index of current postion/ column

    this is support function for calculating transition probabilities.
    this function is only appiled when current column in match state and
    next column is in insertion state.
    """
    trans_counts = transition_dict().copy()
    col = all_columns[ind]
    next_col = all_columns[ind+1]
    col_state = get_states(all_columns, ind)
    next_col_state = get_states(all_columns, ind+1)
            
    while not all_columns[ind][-1] == all_columns[ind+1][-1] =='*':
        ind+=1
        transition = list(map(lambda x,y: x+y, col_state, next_col_state))
        # add to count dict
        for result in transition:
            if result in trans_counts:  # transition happens
                trans_counts[result] += 1

        # update col_state
        for index, state in enumerate(next_col_state):
            if state != '': # if there is insertion
               col_state[index] = state 
        next_col = all_columns[ind+1]
        next_col_state = get_states(all_columns, ind+1)
    return trans_counts

def transition_counts(all_columns):
    """ Return dict of dict, the transition counts of each position/column

    all_columns:list of list, each element is a list of string

    count the number of transitions between states. Return a dict of dict,
    each position as key, and transition counts dict as value. e.g.,
    {0: {'mm':6, 'mi':0,'im':0, 'md':1,..},  1: {'mm':6, 'mi':0...},...} 
    """
    all_counts = {}
    d = 0
    for i in range(len(all_columns)-1):
        if all_columns[i][-1] == '*':
            if all_columns[i+1][-1] == '*': # if next column in match state
                counts = match_match(all_columns, i)
            else: # if next column in insertion state
                counts = match_insertion(all_columns, i)
            all_counts[d] = counts
            d+=1
    return all_counts

def transition_probs(all_counts):
    """return dict of dict of dict, position as key, transition probs dict as
    value
       e.g.,{0: {'m': {'mm': 0.7, 'mi': 0.1, 'md': 0.2},
                'i': {'im': 0.333, 'id': 0.333, 'ii': 0.333},
                'd': {'dm': 0.333, 'dd': 0.333, 'di': 0.333}},...}
                
    all_counts: dict of dict， position as key， dict of transition probs of
    that position as value. e.g.,{0: {'mm':6, 'mi':0,'im':0, 'md':1,..},
                                  1: {'mm':6, 'mi':0,...},...}

    this function calculate the transition probability of each position,
    the probabilities of a certain state to other states should sum up to 1.
    """
    m = ['mm','mi','md']
    i = ['im','id','ii']
    d = ['dm','dd','di']
    trans_probs = {}
    for col_ind in all_counts.keys(): # for each column
        trans_probs[col_ind] = {}
        counts = all_counts[col_ind]
        for trans in counts: # for each transition
            counts[trans] += 1 # add pseudo count 1

    # normalize transition count to probs
        denominator_m = counts['mm']+counts['mi']+counts['md']
        denominator_i = counts['im']+counts['ii']+counts['id']
        denominator_d = counts['dm']+counts['di']+counts['dd']
        dict_m = {}
        for key_m in m:
            dict_m[key_m] = round(counts[key_m]/denominator_m,3)
        trans_probs[col_ind]['m'] = dict_m
        dict_i = {}
        for key_i in i:
            dict_i[key_i] = round(counts[key_i]/denominator_i,3)
        trans_probs[col_ind]['i'] = dict_i
        dict_d = {}
        for key_d in d:
            dict_d[key_d] = round(counts[key_d]/denominator_d,3)
        trans_probs[col_ind]['d'] = dict_d
    return trans_probs       
        
        
def sample(events):
    """Return a key from dict based on the probabilities 

    events: dict of {key: probability}, sum of probabilities
    should be 1.0. 
    """
    k = list(events.keys())
    cum = [0 for i in k]

    cum[0] = events[k[0]]
    for i in range(1,len(events)):
        cum[i] = cum[i-1] + events[k[i]]
    # Should not be necessary, but for safety
    cum[len(cum)-1] = 1.0

    r = random()
    pick = ''
    i = 0
    while (pick == '') and (i < len(cum)):
        if r < cum[i]:
            pick = k[i]
        i = i + 1
    return pick


def seq_generator(emissions, transitions, nrmatch):
    """Return a string, generate random sequence use obtained HMM model

    emissions: dict of dict, int (0,1,..,nrmatch-1) as key, emissions dict
    as value, each key in emissions dict is a string/amino acid, value is a
    float/ probability
    transitions: dict of dict of dict, position as key, transition probs dict as
    value. e.g.,{0: {'m': {'mm': 0.7, 'mi': 0.1, 'md': 0.2},
                     'i': {'im': 0.333, 'id': 0.333, 'ii': 0.333},
                     'd': {'dm': 0.333, 'dd': 0.333, 'di': 0.333}},...}
    nrmatch: int, the total number of match states
    """
    seq = ''
    
    # first transition
    i = 0
    state = 'm'
    seq += sample(emissions[0])
    for i in range(nrmatch-1):
        trans = sample(transitions[i][state])
        state = trans[-1]
        if state == 'm':
            amino_acid = sample(emissions[i+1])
        elif state == 'd':
            amino_acid = ''
        else:
            amino_acid = sample(pa)
        seq += amino_acid
    return seq
        


if __name__ == "__main__":

    # implement main code here
    infile = 'test.fasta'
    data_test = fastaParser(infile)
    file_large = 'test_large.fasta'
    data_large = fastaParser(file_large)

 
    # Question 1. nr of match state
    all_columns = mark_state(data_test)
    columns_large = mark_state(data_large)
    nrmatch = [col[-1] for col in all_columns].count('*')
    nrmatch_large = [col[-1] for col in columns_large].count('*')
    for i in all_columns:
        print(i)
    print('number of match states:', nrmatch)
    

    # Question 2.1 emission probabilities
    emission_probabilities = emission_probs(all_columns)
    emission_large = emission_probs(columns_large)
    print('emission probs for test.fasta')
    for i in range(len(all_columns)):
        print(i,':',emission_probabilities[i])

    # Question 2.2 transition probabilities
    all_counts = transition_counts(all_columns)
    transitions = transition_probs(all_counts)
    print('\ntransition probs for test.fasta')
    for ind in all_counts.keys():
        print(ind,':',transitions[ind])

    counts_large = transition_counts(columns_large)
    transitions_large = transition_probs(counts_large)
    print('\n\ntransition probs for test_large.fasta')
    for ind in counts_large.keys():
        print(ind,':',transitions_large[ind])    

    # Question 3  emission probs in PSSM format
    aa = ['A', 'C', 'D', 'E', 'F', 'G','H', 'I', 'L', 'K', 'M', 'N',\
    'P', 'Q', 'R', 'S', 'T', 'V','W', 'Y']
    print('\nemssion probs for test.fasta in PSSM format\n')
    print(' ', *aa, sep = '      ')
    matrix = [[] for i in range(len(all_columns))]
    for i in range(len(all_columns)):
        for amino_acid in aa:
            matrix[i].append(round(emission_probabilities[i][amino_acid],3))
    for ind, record in enumerate(matrix):
        print(ind, record)
        
    print('\n\nemssion probs for test_large.fasta in PSSM format\n')
    print(' ', *aa, sep = '        ')
    matrix_2 = [[] for i in range(len(columns_large))]
    for i in range(len(columns_large)):
        for amino_acid in aa:
            matrix_2[i].append("%.3f" %emission_large[i][amino_acid])
    for ind, record in enumerate(matrix_2):
        print(ind, record)

    # Question 4
    print('randomly generated sequences for test.fasta\n')
    for i in range(10):
        seq = seq_generator(emission_probabilities, transitions, nrmatch)
        print(i+1,seq)

    print('\n\nrandomly generated sequences for test_large.fasta\n')
    for i in range(10):
        seq = seq_generator(emission_large, transitions_large, nrmatch_large)
        print(i+1,seq)
        
    

    # Question 5
    print('\n\nrandomly generated sequences using background probs\n')
    for i in range(10):
        emission_pa = {}
        for i in range(nrmatch):
            emission_pa[i] = pa
        seq = seq_generator(emission_pa, transitions, nrmatch)
        print(seq)    
    


