#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Xinyuan Min
Student number: 950829573070
Implementation of the SSAHA algorithm

Hints:
- write a function to generate a hash table of k-mers from a set of sequences
- write a function to return a sorted list of hits in the hash table
for a query sequence
- write a function to find the best hit
- write a fasta parser to read in the Arabidopsis data

"""
# import statements
from sys import argv

# implement your functions here

def hashTable(seqs, k):
    """Return a dictionary which contains all k-mers in the subject sequences（key），
    and their corresponding positions (value). The value is a list containing one
    or several tuples.

    seqs: a list of multiple strings (subject sequences)
    k: a positive integer, defines the number of bases in a k-tuple

    this function is to construct a hash table for all k-mers in the subject sequences
    in the database.
    """
    kmer_positions = {}
    # for i, seq in enumerste(seqs):
        # kmer_position[seq[j: j+k].append(i+1, j+1)
    for i in range(len(seqs)):
        for j in range(0, len(seqs[i])-k+1, k):
            # kmer = seqs[i][j:j+k]
            # add the sliced kmer direct to dictionary, this would save the second slicing time
            if seqs[i][j:j+k] not in kmer_positions:
                kmer_positions[seqs[i][j:j+k]] = []
            kmer_positions[seqs[i][j:j+k]].append((i+1, j+1))
    return kmer_positions  


def sortHits(seqs, query, k):
    """Return a list of tuples in qscending order.

    seqs: a list of multiple strings (subject sequences)
    query: string
    k: a positive integer, defines the number of bases in a k-tuple

    to calculate the (index, shift, offset) for each k-mer in the query sequence,
    and sort all the hits in ascending order
    """ 
    hash_table = hashTable(seqs, k)
    hits = []
    for t in range(len(query) -k+1):
        # get the positions from hash table for each k-mer in query sequence 
        if query[t:t+k] in hash_table:
            positions = hash_table[query[t:t+k]]
            # calculate column M in table 2, hits is a list of list
            hit = list(map(lambda x: (x[0], x[-1]- t, x[-1]), positions))
            hits.append(hit)
    # unpack list of list (hits), return sorted hits list
    sorted_hits = []
    for lst in hits:
        for tup in lst:
            sorted_hits.append(tup)
    sorted_hits.sort()
    return sorted_hits

def longestMatch(seqs, query, k):
    """Return a string of the longest match between query and database

    seqs: a list of multiple strings (subject sequences)
    query: string
    k: a positive integer, defines the number of bases in a k-tuple
    
    this function is to get the start and end indexes of matched sequence
    and query sequence, and return the longest match of the two sequence.
    """ 
    hits_list = sortHits(seqs, query, k)
    # get (index, shift) for each hits (without offset)
    index_shift = list(map(lambda x: (x[0], x[1]), hits_list))
    
    counts = []
    # count how many time a specific (index, shift) combination appears
    for tup in index_shift:
        counts.append(index_shift.count(tup))
    # get the index of the most shared (index, shift) combination,
    # e.g. index of (2, 7, 7) in column M in that paper 
    start = counts.index(max(counts))
    match_hits = []
    # get all the hits that shared the same most shared (index, shift)
    # e.g. a list containing 4 highlighted hits in column M in that paper 
    for idx in range(start, len(index_shift)):
        if index_shift[idx] == index_shift[start]:
            match_hits.append(hits_list[idx])
            while index_shift[idx] != index_shift[start]:
                break
    # obtain the matched sequence in the database
    matched_seq = seqs[match_hits[0][0]-1]

    # begin & end are the indexes of start and end base in matched sequence,
    # begin_query & end_query are the indexes of start and end base in query sequence 
    begin = match_hits[0][2] - 1
    end = match_hits[-1][2] -1 + k -1
    begin_query = match_hits[0][2] - match_hits[0][1] 
    end_query = match_hits[-1][2] - match_hits[-1][1] + k -1
    
    # check if the before and next bases of query and matched sequence are identical,
    # if so, add the identitical base to longest match to improve the sensitivity of this algorithm
    if (end_query +1) <= (len(query)-1) and (end+1) <= (len(matched_seq)-1):
        if query[end_query +1] == matched_seq[end+1]:
            end_query += 1
            end += 1            
    if (begin_query -1) >= 0 and (begin-1) >= 0:
        if query[begin_query -1] == matched_seq[begin-1]:
            begin_query -= 1
            begin -= 1

    longest_match = matched_seq[begin:end+1]
    # print the alignment 
    align_up = ' '*(begin-begin_query) + query
    align_down = matched_seq
    print(align_up)
    print(align_down)
    
    return longest_match, match_hits



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






if __name__ == "__main__":
    k = int(argv[1])
    query = 'TGCAACAT'

    
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1, s2, s3]
    
    # Question 1
    #for kmer in sorted(hashTable(seqs, k).keys()):
        #print(kmer, hashTable(seqs, k)[kmer])

    # Question 2
    #sorted_hits = sortHits(seqs, query, k)
    #print('number of hits:', len(sorted_hits), '\n','first hit:', sorted_hits[0], '\n','last hit:', sorted_hits[-1])

    # Question 3
    #print(longestMatch(seqs, query, k))

    # Question 4
    #print(len(fastaParser('TAIR10.fasta.txt')))
    #total_length = 0
    #for seq in fastaParser('TAIR10.fasta.txt'):
        #total_length += len(seq)
    #print('total length:',total_length)

    # Question 5
    #arabidopsis = fastaParser('TAIR10.fasta.txt')
    #print('There are',len(hashTable(arabidopsis,k)), k, '-mers in Arabidopis chromosomes')
    
    # Question 6
    #arabidopsis = fastaParser('TAIR10.fasta.txt')
    #athal_query = fastaParser('athal_query.fasta.txt')
    #aquery = ''.join(athal_query)    
    #print(longestMatch(arabidopsis, aquery,k))
  
