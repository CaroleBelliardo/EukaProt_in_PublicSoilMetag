#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import sys
import multiprocessing as mp
import pandas as pd
import pathlib  # test file
from Bio import SeqIO
from time import process_time
import pyfasta
from modules import gestError
# from gff3 import Gff3

wd = os.path.dirname(os.path.realpath(__file__))
krakenModels_dico=dict()
## -- Parse kraken
def time(func):
    t1_start = process_time()
    func
    t1_stop = process_time()
    print(f'time : {t1_stop - t1_start}')


def progress(iteration, steps, max_value, no_limit=False):
    if int(iteration) == max_value:
        if no_limit == True:
            sys.stdout.write('\r')
            print("[x] \t%d%%" % (100), end='\r')
        else:
            sys.stdout.write('\r')
            print("[x] \t%d%%" % (100))
    elif int(iteration) % steps == 0:
        sys.stdout.write('\r')
        print("[x] \t%d%%" % (float(int(iteration) / int(max_value)) * 100), end='\r')
        sys.stdout.flush()
    else:
        pass


def appendIndicoValue(dico, key, value_list):
    if key in dico:
        dico[key] = dico[key] + value_list
    else:
        dico[key] = value_list
    return dico


def appendSetIndicoValue(dico, key, value_set):
    if key in dico:
        dico[key] = dico[key].union(value_set)
    else:
        dico[key] = value_set
    return dico


def parseLineageIntoList(string):
    liste = list(map(int, (filter(None, string.split(';')))))
    return liste

def parrallelize(func, jobL):
    pool = mp.Pool(os.cpu_count())
    for i, _ in enumerate(pool.imap_unordered(func, jobL), 1):
        progress(i, 1, len(jobL))



def appendDicoIndicoValue(dico, key, key2, value_set):
    if key in dico:
        if key2 in dico[key]:
            dico[key][key2] = dico[key][key2].union(value_set)
        else:
            dico[key][key2] = value_set
    else:
        dico[key] = {key2: value_set}
    return dico


def new_fasta_extract(t_list):
    mid, fastaPath, model_db = t_list
    fastaFile_in_str=gestError.file_exist(f'{fastaPath}/{mid}.a.fna_filtered')
    faa_dico = pyfasta.Fasta(str(fastaFile_in_str))
    for mod in model_db.index  :
        file_p = pathlib.Path(model_db.model_mid_contigs[model_db.index == mod].values[0]+'/'+mid)
        if file_p.exists() :
            df = pd.read_csv(file_p, delimiter=',',names=["seq"])
            seqs=set(df.seq)
            fastaOut=model_db.fnaPath[model_db.index == mod].values[0]+'/'+mid+'.fna'
            with open(fastaOut, 'a') as handle2:
                [handle2.write(f'>{i}_{mid}\n{str(faa_dico[i])}\n') for i in faa_dico if i in seqs]


def extractSeqRun(fastaRepo, Models_df_fasta_contig):  # metagL,ContigK_fasta,taxon,Kraken_taxon,contigs
    list_metag=set()
    for m in Models_df_fasta_contig.model_mid_contigs:
        list_metag = list_metag.union(os.listdir(m))
    jobs=[]
    for m in list_metag :
        jobs.append([m,fastaRepo,Models_df_fasta_contig])
    parrallelize(new_fasta_extract, jobs)


def concatFile(repo):  # cat repo ; rm tmp file
    if len(os.listdir(repo)) != 0:
        str1 = "cat " + repo + "/* > " + repo + ".fna; rm -r " + repo
    else:
        str1 = "rm -r " + repo
    content = os.popen(str1).read()


def checkFile_split(fichier):
    if os.stat(fichier).st_size == 1000000000:
        print('need to be splited')
        # split_fasta fichier > fichier+'_split'

# Augustus output traitement
# analyse gff
def gfftobed(l):  # conversion en fichier bed
    gffP, BedP = l  # ex . : ['test/GFF_krakenEuka-contigs/', 'test/BED_krakenEuka-contigs/', 'test/AA_krakenEuka-contigs/', 'Mammalia.gff']
    GffIN = f'{gffP}.gff'
    BedOUT = f'{BedP}.bed'
    str1 = f'awk -F\'\\t\' \'$3~/^gene/\' {GffIN} | awk -F\'\\t\' ' + "'{print $1,$4,$5,$6,$7,$8,$9,$2,$3}'" + f' OFS=\'\\t\' > {BedOUT}'
    content = os.popen(str1).read()


def gfftofasta(l):  # extraction seq prot
    gffP, aaP = l
    GffIN = f'{gffP}.gff'
    aaOut = f'{gffP}.aa'
    aaMV = f'{aaP}.aa'

    str1 = "perl " + wd + "/getAnnoFasta.pl " + GffIN
    content = os.popen(str1).read()

    if pathlib.Path(aaOut).exists():
        str2 = f'mv {aaOut} {aaMV}'
        content = os.popen(str2).read()

    else:
        print(f'    - error; no "aa" named :{aaOut}')


def gffParse(Gff_Bed_Aa_df,log_f):  # fastaRepo,dico_contigTaxon,fastaOutRepo      #metagL,ContigK_fasta,taxon,Kraken_taxon,contigs
    # listeOfgff = os.listdir(liste[0])
    jobL_Bed = []
    jobL_Fasta = []
    for i in Gff_Bed_Aa_df.index:
        try:
            gff = open(Gff_Bed_Aa_df.gffPath[i]+'.gff' , 'r')
            gff.close()
            jobL_Bed.append([ Gff_Bed_Aa_df.gffPath[i], Gff_Bed_Aa_df.bedPath[i]])#
            jobL_Fasta.append([Gff_Bed_Aa_df.gffPath[i], Gff_Bed_Aa_df.aaPath[i]])#Gff_Bed_Aa_df.gffPath[i]]
        except IOError:
            print(f'    - error; no "gff" named :{Gff_Bed_Aa_df.gffPath[i]}', file=open(log_f, 'a'))
    parrallelize(gfftobed, jobL_Bed)
    parrallelize(gfftofasta, jobL_Fasta)


def addLineage(f):
    str1 = "sh ./modules/addL.sh" + f
    content = os.popen(str1).read()



def filtreSeq(fastaIn, ContigList, fastaOut):
    fasta_sequences = SeqIO.parse(open(fastaIn), 'fasta')  # open fasta
    for seq in fasta_sequences:
        if seq.id in ContigList:  # if sequence is in list of contig of this taxon for this metag
            SeqIO.write([seq], fastaOut, "fasta")  # print sequence in the file

def dmdParse(l):
    aa, dmdout, aaEuk = l
    addLineage(dmdout)
    dmdout = dmdout + '.lineage'
    b = pd.read_table(dmdout, sep='\t', names=['pName', 'taxid', 'evalue', 'lineage'])
    b[['contigid', 'protid']] = b.pName.str.split(".", 1, expand=True, )  # split en 2 col
    b['freqContig'] = b.groupby(by='contigid')['contigid'].transform('count')
    b['freqContige'] = b.groupby(['contigid', 'e'])['contigid'].transform('count')
    bb = b.loc[b['e'] == True]
    l_eukContig = bb[bb.freqContige / bb.freqContig > 0.5].pName.to_list()
    filtreSeq(aa, l_eukContig, aaEuk)
    with open(aaEuk + '.list', 'w') as f:
        for item in l_eukContig:
            f.write("%s\n" % item)

    # extract Prodigal seq
    def extProd(p):
        print('ext prod')
        pl = p + '.list'
        print(pl)
        f = open(pl, "r")
        liste = f.readlines()  # list euka contigs
        for i in liste:
            ii = i.rsplit('_', 1)
            contigId, metagId = ii
            ## ou
        # DT=pd.read_table(args.modelTable, sep='\t');
        ##DK[['contig','metagId']] = DK.contig.str.rsplit("_",1,expand=True,)
