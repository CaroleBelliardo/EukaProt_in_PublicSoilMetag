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

#parser.add_argument('-g', '--gff', type=str, help='GFF directory path with gff files (.gff) \n')
#parser.add_argument('-faa', '--proteins', type=str, help='Proteins directory path with prot files (.faa) \n', default= None)
args = parser.parse_args()


##--- Functions -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def makeH (F): # add header to kraken file
    out = pathlib.Path(F+"_h")
    if not out.exists ():
        str1 = "sed '1s/^/stat\\tcontig\\tf\\ttaxid\\tlist\\tlineage\\n/' "+F+" > "+out # duplique kraken file with header
        content = os.popen(str1).read()
    else :
        print(F+'_h already exit')
        
def mkdir_exist(repoP) :# creation repertoire temporaire
    if not os.path.exists(repoP): 
        os.makedirs(repoP)

def file_exit(filepath,func): # test if file existe
    fileP = pathlib.Path(filepath)
    if fileP.exists ():
        print('!! '+filepath+' ok')
        func(filepath)
    else:
        print (filepath+"file not exist")

# Parse file
def extractSeq(aMetag_fnaPath,uMetag_fnaPath,ContigL,Outputfile):
    afileP = pathlib.Path(aMetag_fnaPath) # test if fna file exist
    ufileP = pathlib.Path(aMetag_fnaPath)
    if afileP.exists ():
        Metag_fnaPath=aMetag_fnaPath
    elif ufileP.exists ():
        Metag_fnaPath=uMetag_fnaPath
    else :
        Metag_fnaPath=None
        print ("error: "+uMetag_fnaPath+' and '+aMetag_fnaPath+" files dont exist")
    if Metag_fnaPath != None :
        fasta_sequences = SeqIO.parse(open(Metag_fnaPath),'fasta') # open fasta
        for seq in fasta_sequences:
            if seq.id in ContigL:    # if sequence is in list of contig of this taxon for this metag 
                with open (Outputfile,'a') as i:
                    SeqIO.write([seq], i, "fasta") # print sequence in the file
        
def extractSeqByMetag(metagid,out): # extraction de sequence
    contigL=Kraken_taxon.contig.loc[Kraken_taxon['metagId'] == metagid].tolist()  # list des contigs du metag == metagID
    ametag_fnaPath=args.contigs+'/'+metagid+'.a.fna' # IN; chemin vers metag
    umetag_fnaPath=args.contigs+'/'+metagid+'.u.fna' # IN; chemin vers metag
    outputfile=out+'/'+metagid+'.fna'    
    extractSeq(ametag_fnaPath,umetag_fnaPath,contigL,outputfile)
    
def krakenToFasta(metagL):
    taxonPath=contigK_fasta+'/'+args.taxon # path variable
    mkdir_exist(taxonPath) # repo specie fasta
# parrallel version
    pool = mp.Pool(mp.cpu_count()) # Step 1: Init multiprocessing.Pool()
    results = [pool.apply(extractSeqByMetag, args=(m,taxonPath)) for m in metagL] # Step 2: `pool.apply` the `howmany_within_range()`   
    pool.close()# Step 3: Don't forget to close


#SING2 SING_IMG 
def augustus(model, metag, gff): ## ajouter singularity avec variable sing et sing2
    str1 = " augustus --uniqueGeneId=true --gff3=on --species="+model+" "+metag+" > "+gff
    content = os.popen(str1).read()
    
def augustusRun(Mid,Model,gffP):
    MetagPath=contigK_fasta+'/'+args.taxon+'/'+Mid+'.fna' 
    ProtPath=gffP+'/'+Mid+'.gff' 
    fileP = pathlib.Path(MetagPath) # test if fna file exist
    if fileP.exists ():
        augustus(Model, MetagPath, ProtPath)
    else:
        print('error: file '+MetagPath+' doesnt exist' )
        
def augustusLoop():
    taxonPathgff=auggff+'/'+args.taxon # path variable
    mkdir_exist(taxonPathgff) # repo specie fasta
#variable for loop
    model=dt.model.loc[dt['deeper_taxon'] == args.taxon].unique()[0]
    lseukfna=os.listdir(contigK_fasta+'/'+args.taxon); lseukfna=[x.strip('.fna') for x in lseukfna]
# parrallel version
    pool = mp.Pool(mp.cpu_count()) # Step 1: Init multiprocessing.Pool()
    results = [pool.apply(augustusRun, args=(m,model,taxonPathgff)) for m in lseukfna] # Step 2: `pool.apply` the `howmany_within_range()`
    pool.close()# Step 3: Don't forget to close
    return(taxonPathgff)

# Augustus output traitement
#analyse gff
def gfftobed(gff,AugProt_bed):# conversion en fichier bed
    str1 = " awk -F'\t' '$3~/^gene/' "+auggff+'/'+gff+".gff | awk -F'\t' '{print $1,$4,$5,$6,$7,$8,$9,$2,$3}' OFS='\t' > "+AugProt_bed+'/'+gff+".bed"
    content = os.popen(str1).read()

def gfftofasta(gff,augProt_faa): # extraction seq prot
    if os.path.exists(auggff+'/'+gff):
        str1 = " getAnnoFasta.pl "+auggff+'/'+gff+".gff ; mv "+auggff+'/'+gff+".gff.aa "+augProt_faa+'/.'
        content = os.popen(str1).read()

def addH(f,i): # add metagid in fasta header
    str1 = "./addH.sh "+f+' '+i
    content = os.popen(str1).read()

def Concat(r):
    str2 = "cat "+r+"/*_hid > "+r+'.faa' ##add name
    content = os.popen(str2).read()
    
def gffParseLoop(): # augustus output traitement
    augProt_bed=args.output+'/BED_AugustusEuka-prot'; mkdir_exist(augProt_bed)
    augProt_bed=augProt_bed+'/'+args.taxon; mkdir_exist(augProt_bed)
    lsAuggff=os.listdir(auggff); lsAuggff=[x.strip('.gff') for x in lsAuggff]
# parrallel version
    pool = mp.Pool(mp.cpu_count()) #BED
    results = [pool.apply(gfftobed, args=(augf,augProt_bed)) for augf in lsAuggff] # rm result
    pool.close()# Step 3: Don't forget to close
    
    augProt_faa=args.output+'/FAA_AugustusEuka-prot'; mkdir_exist(augProt_faa)
    augProt_faa=augProt_faa+'/'+args.taxon; mkdir_exist(augProt_faa)   
    pool = mp.Pool(mp.cpu_count()) # FAA
    results = [pool.apply(gfftofasta, args=(augf,augProt_faa)) for augf in lsAuggff] #rm result
    pool.close()# Step 3: Don't forget to close
    
    ##add name
    lsAugfaa=os.listdir(augProt_faa); lsAugfaa=[x.strip('.aa') for x in lsAugfaa]
    pool = mp.Pool(mp.cpu_count()) # FAA
    results = [pool.apply(addH, args=(augProt_faa,mid)) for mid in lsAugfaa] #rm result
    pool.close()# Step 3: Don't forget to close
    
    return(augProt_faa)







##dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ajout nom metag fasta seq header
# pool fasta par paquet 

#run diamond

def diamond(I,Proteom, dmdOut): ## ajouter singularity avec variable sing et sing2
    str1 = " diamond blastp -p 60 --more-sensitive  -e 0.000001 --outfmt 102 -d /work/cbelliardo/NR_2020_01_diamond.dmnd -q "+ Proteom+"/"+I+".aa -o " + dmdOut+I+".dmd"
    content = os.popen(str1).read()
    
def diamondLoop(Faa,dmdOut):
    taxonPathdmd=dmdOut+'/'+Taxon ; mkdir_exist(taxonPathdmd) # repo specie fasta
    lsAugfaa=os.listdir(Faa); lsAuggff=[x.strip('.aa') for x in lsAugfaa]
# parrallel version
    pool = mp.Pool(mp.cpu_count()) # Step 1: Init multiprocessing.Pool()
    results = [pool.apply(diamond, args=(i,Faa,taxonPathdmd)) for i in lsAugfaa] # Step 2: `pool.apply` the `howmany_within_range()`
    pool.close()# Step 3: Don't forget to close

#appel diamond 
#filtre diamond
#appel busco

##end dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#def makeH (F): # add header to kraken file
#    str1 = "sed '1s/^/stat\\tcontig\\tf\\ttaxid\\tlist\\tlineage\\n/' "+F+">"+F+"_h" # duplique kraken file with header
#    content = os.popen(str1).read()


##--- Files ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Augustus model table 
dt=pd.read_table(args.modelTable, sep='\t');print('step1 : ok import '+args.modelTable) # import model table 
dt=dt.set_index('species') ;print('step2 : ok index '+args.modelTable)# rownames model table 


# Import Kraken output
file_exit(args.kraken,makeH) ;print('step3 : ok header '+args.kraken)
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
mkdir_exist(contigK_fasta) ;print('step6 : ok CreateRepo kraken') # creation repertoire temporaire
krakenToFasta (metagL) ;print('step7 : ok KrakentoFasta') # >> repo fasta

#run fna genome to faa augustus 
auggff=args.output+'/GFF_AugustusEuka-prot'
mkdir_exist(auggff);print('step8 : ok CreateRepo aug')# creation repertoire temporaire
auggff=augustusLoop();print('step9 : ok augustus')
faa=gffParseLoop();print('step10 : ok aug => bed fasta')

## dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#addHConcat(faa) ;print('step11 : ok concat faa files')
#dmdf=args.output+'/DMD_AugustusEuka-prot'; mkdir_exist(auggff);print('step12 : ok CreateRepo dmd')# creation repertoire temporaire
#diamond(faa,dmdf);print('step13 : ok dmd')

##end dev ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# recup 50%avec script puis 3 metriques 

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
#for file in os.listdir():
 #   if os.path.isfile(os.path.join('Mucoromycota', file)):
  #      print(file)
