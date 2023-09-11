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
# use singularity
#SING_IMG="/bighub/hub/people/carole/work-bighub/sing-image/MetagTools.sif "
#SING2="singularity exec --bind /bighub/hub:/bighub/hub:rw --bind /work/cbelliardo:/work/cbelliardo:rw "

# Options
parser = argparse.ArgumentParser(description='Need PYTHON 3 version with packages numpy, pandas, biopython install with :pip install biopython; pip install numpy; pip install pandas; Ex. command line :  python3 Pipeline_euk_v2.py -fna fasta_test -mt tx_sp_mod_nods.txt -k krak.test -t Mucoromycota -o test ')
parser.add_argument('-f', '--fna', type=str, help='Contigs directory path with genomic files (.fna) \n')
parser.add_argument('-mt', '--modelTable', type=str, help='modele - species (.tab file) \n')
parser.add_argument('-k', '--kraken', type=str, help='Kraken output (.tab file) \n')
parser.add_argument('-t', '--taxon', type=str, help='taxon selected for only one taxon (.tab file) \n')
parser.add_argument('-ps', '--printS', type=bool, nargs='?', default=False, help='modele - species (.tab file) \n')
parser.add_argument('-o', '--output', type=str, help='output repository \n',default= 'PipelineOut')
parser.add_argument('-d', '--delet', type=str, help='output repository \n',default=None)

args = parser.parse_args()
#---
# import augusutus_modelsTable 
DT=pd.read_table(args.modelTable, sep='\t');print('step1 : ok import '+args.modelTable) # import model table 
DT=DT.set_index('model') ;print('step2 : ok index '+args.modelTable)# rownames model table 
gestError.mkdir_exist('tempoFiles')


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

## Main -------------------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
#-- print info models -taxid- deeper taxon dispo pour augustus
    if args.printS != False: #print liste espèces modèles si demande utilisateur en input
        pd.set_option('display.max_rows', DT.shape[0]+1)
#-- run pipeline   
    if args.taxon == None : # run only for one model
        print('  No specific taxon, analyse of all taxon*')
#-- List euka contigs
        Dico_contigTaxon=gestError.krakenTodicoTempo(args.kraken,DT);print('step3 : Dico_contigTaxon ok') # paralleliser !!
#-- Paths 
        fastaOutRepo=args.output+'/FNA_krakenEuka-contigs/'; gestError.mkdir_exist(fastaOutRepo) ;print('step4 : ok CreateRepo kraken_fna')
        # creation sous dossier
        dtdt= DT.deeper_taxon.str.replace(' ','_')
        DT['fnaPath'] = fastaOutRepo+dtdt.astype(str); DT.fnaPath.apply(gestError.mkdir_exist);print('step5 : ok CreateSousRepo kraken_fna')
        GffOutRepo=args.output+'/GFF_krakenEuka-contigs/'; gestError.mkdir_exist(GffOutRepo) ;print('step7 : ok CreateRepo kraken_gff')
        DT['gffPath'] = GffOutRepo + dtdt.astype(str) +'.gff'
        BEDOutRepo=args.output+'/BED_krakenEuka-contigs/'; gestError.mkdir_exist(BEDOutRepo) ;print('step9: ok CreateRepo BED')
        DT['bedPath'] = BEDOutRepo + dtdt.astype(str) +'.bed'
        AAOutRepo=args.output+'/AA_krakenEuka-contigs/'; gestError.mkdir_exist(AAOutRepo) ;print('step10 : ok CreateRepo AA')
        DT['aaPath'] = AAOutRepo + dtdt.astype(str) +'.aa'
        DMDOutRepo=args.output+'/DMD_krakenEuka-contigs/'; gestError.mkdir_exist(DMDOutRepo) ;print('step12: ok CreateRepo DMD')
        DT['dmdPath'] = DMDOutRepo + dtdt.astype(str) +'.dmd'
        AaEukPath=args.output+'/AAeuk_krakenEuka-contigs/'; gestError.mkdir_exist(AaEukPath) ;print('step10 : ok CreateRepo AA')
        DT['aaEukPath'] = AaEukPath + dtdt.astype(str) +'.aa'
#-- extract seq euka
        # parse.extractSeqRun(args.fna,Dico_contigTaxon,fastaOutRepo,DT[['fnaPath','deeper_taxon']]) ;print('step6 : ok extract euka seq ')
#-- concatenate fasta: one file by model      
        # parse.parrallelize(parse.concatFile,DT['fnaPath'].tolist()); print('step7 : ok concat fasta files')        
#-- run augustus
        # callTools.augustusLoop2(DT)
#-- parse augfiles    
        # to_lgff=pd.DataFrame(DT, columns= ['gffPath','bedPath', 'AaPath']).values.tolist()
        # parse.parrallelize(parse.gffParse2,to_lgff);print('step11 : ok parse Gff')
 #-- run diamond
        # to_ldmd=pd.DataFrame(DT, columns= ['aaPath', 'dmdPath']).values.tolist()
        # for ilist in to_ldmd : callTools.diamond(ilist)
        # parse.parrallelize(callTools.diamond,to_ldmd); print('step13 : ok diamond')
  #        to_leuk=pd.DataFrame(DT, columns= ['aaPath', 'dmdPath','aaEukPath']).values.tolist()
 #        parse.parrallelize(parse.dmdParse,to_leuk)
# extract prodigal prot
        parse.parrallelize(parse.extProd,DT['aaEukPath'].values.tolist())
        # 

 #        ## add BUSCO!!! 
 #    else : 
        Taxon=args.taxon
        # print(Taxon+' specific taxon analyse*')
        # ... ... ...
