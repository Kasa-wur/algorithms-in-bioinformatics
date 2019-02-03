# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os.path

 
"""
Name: Xinyuan Min, 950829573070
Script of first assignment
"""


file_name = sys.argv[1]
dna_file = open(file_name, "r")
file_content = dna_file.read().replace("\n","")  # get rid of \n at the end of lines
dna_sequence = file_content.upper()  # use upper() to convert all bp to upper case

    
def bpCount(file):
    """Return interger of the length of this DNA sequence

    file: string which represents DNA sequence
    
    This function counts the total amount of base in this DNA sequence
    """
    amount_bp = len(file)
    return amount_bp 
    
def contentGC(file):
    """Return the percentage of GC content with 2 decimals。

    file：string which represents DNA sequence

    This function calculates GC content of this DNA sequence in a percentage
    value with 2 decimals
    """
    g_content = file.count("G")
    c_content = file.count("C")
    content = (g_content + c_content) / bpCount(file)
    # use format() to limit output to 2 decimals
    gcContent = "{:.2%}".format(content)
    return gcContent

def reverseSequence(file):
    """Return string of the reverse sequence

    file：string which represents DNA sequence

    If base in bpdict
    """
    bp_dict = {"A" : "T",
               "T" : "A",
               "G" : "C",
               "C" : "G"
               }
    reverse_sequence = ""
    for base in file:
        complementary_base = bp_dict[base]
        reverse_sequence += complementary_base
    return reverse_sequence

        
if __name__ == "__main__":
    #calculate the length of the sequence
    length = bpCount(dna_sequence)
    print("The length of the sequence is ", length, " bp")
    #calculate the GC content of the sequence
    gccontent = contentGC(dna_sequence) 
    print("The GC content of the sequence is ", gccontent)
    # return the reverse complemant of the sequence
    reverse_sequence = reverseSequence(dna_sequence)
    print("This is the reverse complement of the sequence:\n",reverse_sequence)

          
