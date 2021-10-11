#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import sys
#from modules import gestError
#from modules import callTools
import multiprocessing as mp
#import pandas as pd
#import pathlib # test file
from Bio import SeqIO
#import operator
import pyfastx

def parrallelize(func,jobL):   
    pool=mp.Pool(60)
    for i, _ in enumerate(pool.imap_unordered(func,jobL),1):
        progress(i,1,len(jobL))  

        

def extractSeq(t_list):
    seqId_set, fastaPath,outputFasta= t_list
    fasta_sequences = SeqIO.parse(open(fastaPath),'fasta') # open fasta
    output_handle = open(outputFasta, "a")
    for seq in fasta_sequences:
        if seq.id in seqId_set :
            SeqIO.write([seq], output_handle , "fasta") # print sequence in the file
    output_handle.close()
    

    
