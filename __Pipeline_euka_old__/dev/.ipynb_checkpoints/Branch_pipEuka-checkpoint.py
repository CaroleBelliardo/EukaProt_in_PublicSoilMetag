#!/usr/bin/python3
# -*- coding: utf-8 -*-
#TEST
import pprint   # debugg -- affichage Dump des structures de données
import os       # manip systeme
import pathlib # test file
import sys      # gestion arguments et opt
import argparse
import numpy as np
import pandas as pd
import re
from modules import gestError
from modules import parse
from modules import callTools
from time import process_time

t1_start = process_time() 

# Options
parser = argparse.ArgumentParser(description='Need PYTHON 3 version with packages numpy, pandas, biopython install with :pip install biopython; pip install numpy; pip install pandas; Ex. command line :  python3 Pipeline_euk_v2.py -fna fasta_test -mt tx_sp_mod_nods.txt -k krak.test -t Mucoromycota -o test ')
parser.add_argument('-f', '--fna', type=str, help='Contigs directory path with genomic files (.fna) \n')
parser.add_argument('-mt', '--modelTable', type=str, help='modele - species (.tab file) \n')
parser.add_argument('-k', '--kraken', type=str, help='Kraken output (.tab file) \n')
parser.add_argument('-l', '--lineage', type=str, help='lineage file, (.tab file) \n',default= '/bighub/hub/DB/ncbi_taxo/20200724/20200724_taxonomy_lineage.txt')
parser.add_argument('-o', '--output', type=str, help='output repository \n',default= 'PipelineOut')
parser.add_argument('-d', '--delet', type=str, help='output repository \n',default=None)
parser.add_argument('-e', '--eukTaxid', type=str, help='based taxid looking for \n',default=2759)
parser.add_argument('-m', '--defaultModel', type=str, help='based taxid looking for \n',default="leishmania_tarentolae")

args = parser.parse_args()
#---
gestError.mkdir_exist('tempoFiles')


## Main -------------------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
## Kraken to fna datasets
#-- import augusutus_modelsTable
    gestError.file_exist(args.modelTable)
    Models_df=pd.read_table(args.modelTable, sep='\t');print('step1 : ok import '+args.modelTable) # import model table 
    Models_df=Models_df.set_index('model') ;print('step2 : ok index '+args.modelTable)# rownames model table
   
#-- parse kraken output
    gestError.file_exist(args.kraken)
    KkOutput_df=pd.read_table(args.kraken, sep='\t', names=['classe','seqid','taxid','seqLen','LCAListe']);print('step3 : ok import '+args.kraken) # import model table
    KkOutput_df=parse.unclassFilter(KkOutput_df)    # unclassified seq
    KkOutput_df=gestError.nonEuk(KkOutput_df) # already non euka listed seq

#-- List euka contigs
    Dico_contigTaxon=gestError.krakenTodicoTempo(Models_df,KkOutput_df,args.lineage,args.eukTaxid,args.defaultModel)
    # => {model: ['Ga0272421_1202915_3300031909'], 10090: ['Ga0272421_1213909_3300031909'] ...
    # print(Dico_contigTaxon)
#-- Paths 
    fastaOutRepo=args.output+'/FNA_krakenEuka-contigs/'; gestError.mkdir_exist(fastaOutRepo) ;print('step8 : ok CreateRepo kraken_fna')
    # creation sous dossier
    dtdt= Models_df.deeper_taxon.str.replace(' ','_')
    Models_df['fnaPath'] = fastaOutRepo+dtdt.astype(str); Models_df.fnaPath.apply(gestError.mkdir_exist);print('step9 : ok CreateSousRepo kraken_fna')      
    GffOutRepo=args.output+'/GFF_krakenEuka-contigs/'; gestError.mkdir_exist(GffOutRepo) ;print('step7 : ok CreateRepo kraken_gff')
    Models_df['gffPath'] = GffOutRepo + dtdt.astype(str) +'.gff'
    BEDOutRepo=args.output+'/BED_krakenEuka-contigs/'; gestError.mkdir_exist(BEDOutRepo) ;print('step9: ok CreateRepo BED')
    Models_df['bedPath'] = BEDOutRepo + dtdt.astype(str) +'.bed'
    AAOutRepo=args.output+'/AA_krakenEuka-contigs/'; gestError.mkdir_exist(AAOutRepo) ;print('step10 : ok CreateRepo AA')
    Models_df['aaPath'] = AAOutRepo + dtdt.astype(str) +'.aa'
    DMDOutRepo=args.output+'/DMD_krakenEuka-contigs/'; gestError.mkdir_exist(DMDOutRepo) ;print('step12: ok CreateRepo DMD')
    Models_df['dmdPath'] = DMDOutRepo + dtdt.astype(str) +'.dmd'
    AaEukPath=args.output+'/AAeuk_krakenEuka-contigs/'; gestError.mkdir_exist(AaEukPath) ;print('step10 : ok CreateRepo AA')
    Models_df['aaEukPath'] = AaEukPath + dtdt.astype(str) +'.aa'
# #-- extract seq euka
    # parse.extractSeqRun(args.fna,Dico_contigTaxon,fastaOutRepo,Models_df[['fnaPath','deeper_taxon']]) ;print('step10 : ok extract euka seq ')
# #-- concatenate fasta: one file by model      
    parse.parrallelize(parse.concatFile,Models_df['fnaPath'].tolist()); print('step11 : ok concat fasta files')
# # -- dev ----***
# #-- check file size
#     ## MODIF. => Pas de creation repo mais split file dans Repo courantet ajout a la liste des fichier;
#     # Pour chaque modekl, PAth variable == List!!
#     # ajout iteration sur liste repo à toutes les etapes !! 
#         # parse.parrallelize(parse.checkFile_split,Models_df['fnaPath'].tolist()); print('step7.5 : ok check fasta files size')
# # -- dev end ----***        
#-- run augustus
    callTools.augustusLoop2(Models_df)
# -- parse augfiles    
    to_lgff=pd.DataFrame(Models_df, columns= ['gffPath','bedPath', 'aaPath']).values.tolist()
    to_lgff=[GffOutRepo,BEDOutRepo,AAOutRepo]
    parse.gffParse2(to_lgff);print('step11 : ok parse Gff')
 #-- run diamond
    to_ldmd=pd.DataFrame(Models_df, columns= ['aaPath', 'dmdPath']).values.tolist()
    for ilist in to_ldmd : callTools.diamond(ilist)  
    to_leuk=pd.DataFrame(Models_df, columns= ['aaPath', 'dmdPath','aaEukPath']).values.tolist()
    parse.parrallelize(parse.dmdParse,to_leuk)
# # extract prodigal prot.
#         # print(parse.extProd,Models_df['aaEukPath'])
#         # print(parse.extProd,Models_df['aaEukPath'].values.tolist())
#         # parse.parrallelize(parse.extProd,Models_df['aaEukPath'].values.tolist())
#         # 
# 
t1_stop = process_time() 
print("\n\n---\nElapsed time:", t1_start, '->',t1_stop)  
print("Elapsed time during the whole program in seconds:", 
                                         t1_stop-t1_start) 