#!/usr/bin/python3
# -*- coding: utf-8 -*-

##--- Header------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Imports Packages
import pprint   # debugg -- affichage Dump des structures de données
import os       # manip systeme
import pathlib # test file
import sys      # gestion arguments et opt
import argparse
import numpy as np
import pandas as pd
import re
from Bio import SeqIO
import multiprocessing as mp

# Options
parser = argparse.ArgumentParser(description='Need PYTHON 3 version with packages numpy, pandas, biopython,multiprocessing, os, sys, argse, re install with :pip install biopython; pip install numpy; pip install pandas; Ex. command line :  python3 Pipeline_euk.py -fna fasta_test -mt tx_sp_mod_nods.txt -k krak.test -t Mucoromycota -o test ')
parser.add_argument('-fna', '--contigs', type=str, help='Contigs directory path with genomic files (.fna) \n')
parser.add_argument('-g', '--gff', type=str, help='GFF directory path with gff files (.gff) \n')
parser.add_argument('-faa', '--proteins', type=str, help='Proteins directory path with prot files (.faa) \n', default= None)
parser.add_argument('-mt', '--modelTable', type=str, help='modele - species (.tab file) \n')
parser.add_argument('-k', '--kraken', type=str, help='Kraken output (.tab file) \n')
parser.add_argument('-t', '--taxon', type=str, help='taxon selected for only one taxon (.tab file) \n')
parser.add_argument('-ps', '--printS', type=bool, nargs='?', default=False, help='modele - species (.tab file) \n')
parser.add_argument('-o', '--output', type=str, help='output repository \n',default= 'PipelineOut')
args = parser.parse_args()

# variables singularity
SING_IMG="/bighub/hub/people/carole/work-bighub/sing-image/MetagTools.sif"
SING2="singularity exec --bind /bighub/hub:/bighub/hub:rw --bind /work/cbelliardo:/work/cbelliardo:rw"
##--- Functions -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def mkdir_exist(repoP) :# creation repertoire temporaire
    if not os.path.exists(repoP): 
        os.makedirs(repoP)

def makeH (F): # add header to kraken file
    str1 = "sed '1s/^/stat\\tcontig\\tf\\ttaxid\\tlist\\tlineage\\n/' "+F+">"+F+"_h" # duplique kraken file with header
    content = os.popen(str1).read()

def file_exit(filepath,func): # test if file existe
    fileP = pathlib.Path(filepath)
    if fileP.exists ():
        print(filepath+' ok')
        func(filepath)
    else:
        print (filepath+"file not exist")


# Parse file
def extractSeq(Metag_fnaPath,ContigL,Outputfile):
    fasta_sequences = SeqIO.parse(open(Metag_fnaPath),'fasta') # open fasta
    for seq in fasta_sequences:
        if seq.id in ContigL:    # if sequence is in list of contig of this taxon for this metag 
            with open (Outputfile,'a') as i:
                SeqIO.write([seq], i, "fasta") # print sequence in the file

def extractSeqByMetag(kraken_taxon,taxon,metagid,out): # extraction de sequence
    contigL=kraken_taxon.contig.loc[kraken_taxon['metagId'] == metagid].tolist()  # list des contigs du metag == metagID
    metag_fnaPath=args.contigs+'/'+metagid # IN; chemin vers metag
    outputfile=out+'/'+metagid    
    extractSeq(metag_fnaPath,contigL,outputfile)


def krakenToFasta(Taxon):
    taxonPath=contigK_fasta+'/'+Taxon # path variable
    mkdir_exist(taxonPath) # repo specie fasta
    Kraken_taxon=(dk[dk.lineage.str.contains(Taxon)]) # select sous table avec taxon #sst.contig.to_csv(r'tmp/kraken_split/krakenContig_'+Taxon+'.txt',  sep='\t', mode='a', header = None,index =None    
    metagL=Kraken_taxon['metagId'].unique()# list metag with unique occurence
# parrallel version
    pool = mp.Pool(mp.cpu_count()) # Step 1: Init multiprocessing.Pool()
    results = [pool.apply(extractSeqByMetag, args=(Kraken_taxon,Taxon,m,taxonPath)) for m in metagL] # Step 2: `pool.apply` the `howmany_within_range()`   
    pool.close()# Step 3: Don't forget to close
#>> Fasta

##dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#augustus

#SING2 SING_IMG 
def augustus(model, metag, proteome): ## ajouter singularity avec variable sing et sing2
    #sing
    #sing2
    str1 = "augustus --progress=true  --gff3=on --species="+model+" "+metag+" > "+proteome
    content = os.popen(str1).read()

def augustusRun(Mid,Model,Taxon,faaP):
    MetagPath=contigK_fasta+'/'+Taxon+'/'+Mid
    ProtPath=faaP+'/'+Mid
    augustus(Model, MetagPath, ProtPath)
    
def augProtpred(Taxon):
    taxonPathfaa=args.proteins+'/'+Taxon # test if repo prot aug exists for this specie
    mkdir_exist(taxonPathfaa) # repo specie fasta
    #variable for loop
    Kraken_taxon=(dk[dk.lineage.str.contains(Taxon)]) # select sous table avec taxon #sst.contig.to_csv(r'tmp/kraken_split/krakenContig_'+Taxon+'.txt',  sep='\t', mode='a', header = None,index =None    
    metagL=Kraken_taxon['metagId'].unique()# list metag with unique occurence
    model="maize"
    print(dt)
# parrallel version
    pool = mp.Pool(mp.cpu_count()) # Step 1: Init multiprocessing.Pool()
    results = [pool.apply(augustusRun, args=(m,model,taxon,taxonPathfaa)) for m in metagL] # Step 2: `pool.apply` the `howmany_within_range()`
    pool.close()# Step 3: Don't forget to close

augProtpred(args.taxon)
#awk -F'\t' '$3~/^gene/' ${model}_${PREF}.aug > ${model}_${PREF}.gff_gene

#def makeH (F): # add header to kraken file
#    str1 = "sed '1s/^/stat\\tcontig\\tf\\ttaxid\\tlist\\tlineage\\n/' "+F+">"+F+"_h" # duplique kraken file with header
#    content = os.popen(str1).read()

##end dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##--- Files ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mkdir_exist('PipelineOut')

# Augustus model table 
dt=pd.read_table(args.modelTable, sep='\t') # import model table 
dt = dt.set_index('species') # rownames model table

# Import Kraken output
file_exit(args.kraken,makeH) 
dk = pd.read_table(args.kraken+"_h", skipinitialspace=True, usecols=['contig', 'lineage'], sep='\t') # import kraken file
dk[['contig','metagId']] = dk.contig.str.rsplit("_",1,expand=True,) #split contig /metagId

# Kraken to fasta
contigK_fasta=args.output+'/FNA_krakenEuka-contigs'
if (args.proteins == None) : args.proteins=args.output+'/FAA_AugustusEuka-prot'
mkdir_exist(contigK_fasta)# creation repertoire temporaire
krakenToFasta (args.taxon) # >> repo fasta




## Main -------------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    if args.printS != False: #print liste espèces modèles
        pd.set_option('display.max_rows', dt.shape[0]+1)
        print(dt)
    if args.taxon == None :
        print('No specific taxon, analyse of all taxon*')
        # dérouler boucle sur tous les taxons
        # for t in dt.species: 
            #do all analyses Pipeline (t)
    else :
        sst=dt.loc[dt['deeper_taxon'] == args.taxon]
        if sst.empty:
            print(args.taxon+' Not in taxon list, use -ps postion to see all available taxon_node !!!')        
        else:
            print(args.taxon+' in list of taxon***')    
            #do all analyses Pipeline (t)
    
