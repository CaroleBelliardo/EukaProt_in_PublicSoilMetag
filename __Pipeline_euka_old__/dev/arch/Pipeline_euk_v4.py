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
import modules
from modules import gestError
from modules import parse
from modules import callTools

# use singularity
#SING_IMG="/bighub/hub/people/carole/work-bighub/sing-image/MetagTools.sif "
#SING2="singularity exec --bind /bighub/hub:/bighub/hub:rw --bind /work/cbelliardo:/work/cbelliardo:rw "


# Options
parser = argparse.ArgumentParser(description='Need PYTHON 3 version with packages numpy, pandas, biopython install with :pip install biopython; pip install numpy; pip install pandas; Ex. command line :  python3 Pipeline_euk_v2.py -fna fasta_test -mt tx_sp_mod_nods.txt -k krak.test -t Mucoromycota -o test ')
parser.add_argument('-fna', '--contigs', type=str, help='Contigs directory path with genomic files (.fna) \n')
parser.add_argument('-mt', '--modelTable', type=str, help='modele - species (.tab file) \n')
parser.add_argument('-k', '--kraken', type=str, help='Kraken output (.tab file) \n')
parser.add_argument('-t', '--taxon', type=str, help='taxon selected for only one taxon (.tab file) \n')
parser.add_argument('-ps', '--printS', type=bool, nargs='?', default=False, help='modele - species (.tab file) \n')
parser.add_argument('-o', '--output', type=str, help='output repository \n',default= 'PipelineOut')
parser.add_argument('-d', '--delet', type=str, help='output repository \n',default=None)

args = parser.parse_args()


        
##dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ajout nom metag fasta seq header
# pool fasta par paquet 
#run diamond
#filtre diamond
#appel busco

#addHConcat(faa) ;print('step11 : ok concat faa files')
#dmdf=args.output+'/DMD_AugustusEuka-prot'; mkdir_exist(auggff);print('step12 : ok CreateRepo dmd')# creation repertoire temporaire
#diamond(faa,dmdf);print('step13 : ok dmd')

##end dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# recup 50%avec script puis 3 metriques 
# Augustus model table 
dt=pd.read_table(args.modelTable, sep='\t');print('step1 : ok import '+args.modelTable) # import model table 
dt=dt.set_index('species') ;print('step2 : ok index '+args.modelTable)# rownames model table 

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
            gestError.file_exit(args.kraken,gestError.makeH) ;print('step3 : ok header '+args.kraken)
            dk=pd.read_table(args.kraken+"_h", skipinitialspace=True, usecols=['contig', 'lineage'], sep='\t') ;print('step4 : ok import '+args.kraken) # import kraken file
            dk[['contig','metagId']] = dk.contig.str.rsplit("_",1,expand=True,);print('step5 : ok split contig/metagId '+args.kraken) #split contig /metagId
            Kraken_taxon=(dk[dk.lineage.str.contains(args.taxon, na=False)]) # select sous table avec taxon #sst.contig.to_csv(r'tmp/kraken_split/krakenContig_'+Taxon+'.txt',  sep='\t', mode='a', header = None,index =None    
            if args.delet != None :
               d = pd.read_table(args.delet) # tobeContinued
               d=d.mid.astype(str)
               Kraken_taxon=Kraken_taxon[~Kraken_taxon.metagId.isin(d)]
            
            del dk
            metagl=Kraken_taxon['metagId'].unique()# list metag with unique occurence
            metagL=metagl[metagl != np.array(None)]
            
            
            # Kraken to fna fasta 
            contigK_fasta=args.output+'/FNA_krakenEuka-contigs'
            gestError.mkdir_exist(contigK_fasta) ;print('step6 : ok CreateRepo kraken') # creation repertoire temporaire
            parse.krakenToFasta (metagL,contigK_fasta,args.taxon,Kraken_taxon,args.contigs) ;print('step7 : ok KrakentoFasta') #  >>repo fasta
            
            #run fna genome to faa augustus 
            auggff=args.output+'/GFF_AugustusEuka-prot/'+args.taxon
            gestError.mkdir_exist(auggff);print('step8 : ok CreateRepo aug')# creation repertoire temporaire
            auggff=callTools.augustusLoop(auggff,args.taxon,dt,contigK_fasta);print('step9 : ok augustus')
            faa=callTools.gffParseLoop(auggff,args.output,args.taxon);print('step10 : ok aug => bed fasta')
            
      