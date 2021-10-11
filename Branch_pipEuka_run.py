#!/usr/bin/python3
# -*- coding: utf-8 -*-

import pprint  # debugg -- affichage Dump des structures de données
import os  # manip systeme
import argparse
import pathlib
import pandas as pd
from modules import gestError
from modules import parse
from modules import callTools
from time import process_time
import multiprocessing as mp
#from pprint import pprint
#import psutil
import sys

# print(f'{psutil.virtual_memory().percent=}')
# print(f'{sys.getsizeof(st)}')

# Options
parser = argparse.ArgumentParser(
    description='Need PYTHON 3 version with packages numpy, pandas, biopython install with :pip install biopython; '
                'pip install numpy; pip install pandas; Ex. command line :   python Branch_pipEuka.py -f files_test/3300031909.a.fna_filtered -k files_test/krak.txt -t tx_sp_mod_nods.txt -l files_test/taxonomy_lineage.txt_head -o testDev ')
parser.add_argument('-d', '--delet', type=str, help='output repository \n', default=None)
parser.add_argument('-e', '--eukTaxid', type=str, help='based taxid looking for \n', default=2759)
parser.add_argument('-t', '--modelTable', type=str, help='modele - species (.tab file) \n',
                    default='/lerins/hub/DB/Metagenomics/tx_sp_mod_nods.txt')
parser.add_argument('-l', '--lineage', type=str, help='lineage file, (.tab file) \n',
                    default='/lerins/hub/DB/TAXONOMY/online/taxonomy_lineage.txt')

parser.add_argument('-k', '--kraken', type=str, help='Kraken output (.tab file) \n')
parser.add_argument('-f', '--fna', type=str, help='Contigs directory path with genomic files (.fna) \n')
parser.add_argument('-o', '--output', type=str, help='output repository \n', default='EukPip')

args = parser.parse_args()
# ---
gestError.mkdir_exist(args.output)
tmp_r = f'{args.output}/tempoFiles'
gestError.mkdir_exist(tmp_r)
log_f = 'eukPip.log'

## func
def progress(iteration, steps, max_value, no_limit=False):
    global log_f
    if int(iteration) == max_value:
        if no_limit == True:
            sys.stdout.write('\r')
            print("[x] \t%d%%" % (100), end='\r', file=open(log_f, 'a'))
        else:
            sys.stdout.write('\r')
            print("[x] \t%d%%" % (100), file=open(log_f, 'a'))
    elif int(iteration) % steps == 0:
        sys.stdout.write('\r')
        print("[x] \t%d%%" % (float(int(iteration) / int(max_value)) * 100), end='\r', file=open(log_f, 'a'))
        sys.stdout.flush()
    else:
        pass


def parrallelize(func, jobL):
    pool = mp.Pool(os.cpu_count())
    for i, _ in enumerate(pool.imap_unordered(func, jobL), 1):
        progress(i, 1, len(jobL))


def step2(kraken_subt):  # step2
    global Models_df
    mid = kraken_subt.metag.iloc[0]
    for n in Models_df.index.to_list():
        tx = Models_df.taxid_new_deeperTaxon[Models_df.index == n].values[0]
        pathOut = Models_df.model_mid_contigs[Models_df.taxid_new_deeperTaxon == tx].values[0] + '/' + mid
        contigs_df = kraken_subt.seq[kraken_subt.nodes.str.contains(str(tx))]
        if len(contigs_df) > 0: contigs_df.to_csv(f'{pathOut}', index=False, header=False, sep='\t', mode='a')
        kraken_subt = kraken_subt[~(kraken_subt.seq.isin(contigs_df))]  # MAJ kraken
    return(kraken_subt)


def step4(subTab, Nodes_model_dico, metag_id):
    global Ordered_model_byWeight
    global Models_df
    left_taxid_sunT = subTab.taxid.unique().tolist()
    for tx_kraken in left_taxid_sunT:
        kraken_compl_nodes = subTab.nodes[subTab.taxid == tx_kraken].values[0][13:-1].split(';')
        inter = {v: set(k[13:-1].split(';')).intersection(kraken_compl_nodes) for k, v in Nodes_model_dico.items() if
                 len(set(k[13:-1].split(';')).intersection(kraken_compl_nodes)) > 0}
        if len(inter) > 0:
            m = list({k: inter[k] for k in sorted(inter, key=lambda k: len(inter[k]), reverse=True)}.keys())[0]
            if len(m) > 0:
                pathOut = Models_df.model_mid_contigs[Models_df.index == m].values[0] + '/' + metag_id
                newC_list = subTab.seq[subTab.taxid == tx_kraken]
                newC_list.to_csv(f'{pathOut}', index=False, header=False, sep='\t', mode='a')

def step1(kraken_subt):
    t1_start = process_time()
    global Models_df
    global deeperTx_list
    global log_f
    mid = kraken_subt.metag.iloc[0]
    taxid_kraken = set(kraken_subt.taxid)
    intersect_KM = taxid_kraken.intersection(deeperTx_list)
    for tx in intersect_KM:
        pathOut = Models_df.model_mid_contigs[Models_df.taxid_new_deeperTaxon == tx].values[0] + '/' + mid
        contigs_df = kraken_subt.seq[kraken_subt.taxid == tx]
        contigs_df.to_csv(f'{pathOut}', index=False, header=False, sep='\t')
    kraken_subt = kraken_subt[~(kraken_subt.taxid.isin(deeperTx_list))]  # MAJ kraken
    t1_stop = process_time()
    print(f'step 1/3 mid:_{mid}_  krakenToDico :  ok!  taxid == deepertaxid\n     time :   {(t1_stop-t1_start)=}\n\nrunning step 2...', file=open(log_f+'.tmp', 'a'))

#step2
    t1_start = process_time()
    kraken_subt = step2(kraken_subt)
    t1_stop = process_time()
    print(f'step 2/3 mid:_{mid}_ krakenToDico :  run  ok\n     time :   {(t1_stop-t1_start)=}\n\nrunning step 3...', file=open(log_f+'.tmp', 'a'))

#step 3
    t1_start = process_time()

    model_byWeight = dict()
    for model_mid_contigs in Models_df.model_mid_contigs :
        ld=os.listdir(model_mid_contigs)
        model=os.path.basename(model_mid_contigs)
        if len(ld) > 0:
            ln_ld= [len(open(model_mid_contigs+'/'+f).readlines())for f in ld ]
            model_byWeight[model]=sum(ln_ld)
        else:
            model_byWeight[model]=0

    Ordered_model_byWeight=sorted(model_byWeight, key=lambda k: model_byWeight[k], reverse=True)
    nodes_model_dico = {Models_df.complet_taxid[Models_df.index == i].values[0]: i for i in Ordered_model_byWeight}
    step4(kraken_subt,nodes_model_dico,mid)
    t1_stop = process_time()
    print(f'step 3/3 mid:_{mid}_ krakenToDico :  run  ok\n     time :   {(t1_stop-t1_start)=}\n\n\n', file=open(log_f+'.tmp', 'a'))

#
# def appendDicoIndicoValue(dico, key, key2, v_str):
#     if key in dico.keys():
#         if key2 in dico[key]:
#             dico[key][key2]=dico[key][key2]+[v_str]
#             t=dico[key][key2]
#         else:
#             dico[key][key2] = [v_str]
#     else:
#         dico[key] = {key2: [v_str]}
#     return dico
#
# def format_kraken_dico_Fromlist(kraken_model_list):
#     km_dico_formated = dict()
#     for model, contig_list in kraken_model_list:
#         if len(contig_list) > 0 :
#             for c in contig_list :
#                 m=c.rsplit('_',1)[1]
#                 contig=c.rsplit('_',1)[0]
#                 km_dico_formated = appendDicoIndicoValue(km_dico_formated, m, model, contig)
#     return km_dico_formated  # ex:    km_dico={'fileB':{'rhizopus_oryzae': {'Ga0272421_100ustr', 'Ga0272421_100usta', 'Ga0272421_1000028', 'Ga0272421_100ustc', 'Ga0272421_1000035', 'Ga0272421_1000029'}}}
#


# # Main
# -------------------------------------------------------------------------------------------------------------------------------------------------------------
print('Strating ...', file = open(log_f, 'w'))
T_start = process_time()
if __name__ == "__main__":
    krakenModels_dico = dict()
    print('Strating main...', file = open(log_f, 'a'))

    # -- Kraken to fna datasets
#    gestError.file_exist(args.modelTable)  # -- import augusutus_modelsTable
    Models_df = pd.read_table(args.modelTable, sep='\t').set_index('model')

#    print(f'step 1/9 : ok index  {args.modelTable}\n\n', file = open(log_f, 'a'))  # rownames model table

    # -- List euka contigs
#    print('Strating kraken Open...', file=open(log_f, 'a'))

#    km_f = tmp_r + "/krakenModels.txt"
    km_r = tmp_r + "/krakenModels/"
    if pathlib.Path(km_r).exists():  # skip kraken parse
        print('  *** ' + km_r + ' already exist')
    else:
        gestError.mkdir_exist(km_r)
        Models_df['model_mid_contigs'] = [km_r + i for i in Models_df.index]
        Models_df.model_mid_contigs.apply(gestError.mkdir_exist)
        print(f'  *** {km_r} created', file=open(log_f, 'a'))

        kraken_lineage_df = gestError.KrakenLineage_exist(args.kraken, args.lineage, tmp_r, log_f)
        print(f'main : KRAKEN_LINEAGE Imported ok !!! ', file=open(log_f, 'a'))
        # -- all contigsNames with taxid == only 'eukaryote' => print in file; no model could be assign for augustus
        print(f'main : euka taxid seq in file  ... ', file=open(log_f, 'a'))
        euk = kraken_lineage_df[kraken_lineage_df.taxid == args.eukTaxid]
        euk.to_csv(f'{tmp_r}/euk_{args.eukTaxid}.csv', index=False, sep='\t')
        print(f'main : euka taxid seq in file  printed ok  ', file=open(log_f, 'a'))
        kraken_lineage_df = kraken_lineage_df[~(kraken_lineage_df.taxid == args.eukTaxid)]  # remove contigsNames == eukaryote
        print(f'main : euka taxid seq remove from df ok  ', file=open(log_f, 'a'))
        # -- select left row contain euka taxid
        kraken_lineage_df = kraken_lineage_df[kraken_lineage_df['nodes'].notna()] ## A DEPPLACER !!!
        kraken_lineage_df = kraken_lineage_df[kraken_lineage_df['nodes'].str.contains(str(args.eukTaxid))]

        list_of_metag = kraken_lineage_df.metag.unique()
        subtable = [kraken_lineage_df[kraken_lineage_df.metag == t] for t in list_of_metag]
        del kraken_lineage_df

        #  krakenModels_dico = krakenToDico(Models_df, kraken_lineage_df)
        print(f'\n\n\nmain :KRAKEN TO DICO strating ... ', file=open(log_f, 'a'))

        # step 1:  taxid given by kraken == taxid model species; exact search
        t1_start = process_time()
        deeperTx_list = set(Models_df.taxid_new_deeperTaxon)
        parrallelize(step1, subtable)
        del (deeperTx_list)
        t1_stop = process_time()
        print(f'step 2/3: Dico to Kraken ok\n     time :   {(t1_stop-t1_start)=}', file=open(log_f, 'a'))

    repo = {
        'model_mid_contigs' : tmp_r + "/krakenModels/",
        'fnaPath': args.output + '/FNA_krakenEuka-contigs/',
        'gffPath': args.output + '/GFF_krakenEuka-contigs/',
        'bedPath': args.output + '/BED_krakenEuka-contigs/',
        'aaPath': args.output + '/AA_krakenEuka-contigs/',
        'dmdPath': args.output + '/DMD_krakenEuka-contigs/',
        'aaEukPath': args.output + '/AAeuk_krakenEuka-contigs/'}


#    for r in repo.values(): gestError.mkdir_exist(r)
#    print(f'step 3/9 : ok  output Repository created {repo["fnaPath"]} \n\n', file=open(log_f, 'a'))
    # mkdir SUB dossier
    for k, r in repo.items(): Models_df[k] = [r + i for i in Models_df.index]
#    Models_df.fnaPath.apply(gestError.mkdir_exist)
#    print(f'step 4/9 : ok kraken_fna sub Repo {len(repo["fnaPath"])} \n\nrunning extract seq....', file=open(log_f, 'a'))
    #-- extract seq euka
    # list les metag set fichier dans repo model => pool donne metag, Model,
    # appel fonction pour chaque model print fasta si fichier metag exist


#    parse.extractSeqRun(args.fna,  Models_df[['fnaPath','model_mid_contigs']])#
#    print(f'step 5/9 :  ok extract euka seq \n\n', file=open(log_f, 'a'))
# -- concatenate fasta: one file by model
#    parse.parrallelize(parse.concatFile,Models_df['fnaPath'].tolist());
#    print(f'step 6/9 :  ok concat fasta files \n\n', file=open(log_f, 'a'))
#
# #-- run augustus
#    callTools.augustusLoop2(Models_df[['fnaPath','gffPath']])

# -- parse augfiles
    gff_bed_aa_df=Models_df[['gffPath','bedPath', 'aaPath']]
    parse.gffParse(gff_bed_aa_df,log_f)
    print(f'step 7/9 :  ok parse Gff \n\n', file=open(log_f, 'a'))
#
#  #-- run diamond
    list_aa=[i.split('.')[0] for i in os.listdir(repo['aaPath'])]
    to_ldmd=[[Models_df.aaPath[i],Models_df.dmdPath[i]] for i in list_aa ]
    for ilist in to_ldmd : callTools.diamond(ilist)
    del to_ldmd
    t1_start = process_time()
    print(f'step 8/9 :  ok diamond run \n\n', file=open(log_f, 'a'))
    to_leuk = [[Models_df.aaPath[i], Models_df.dmdPath[i], Models_df.aaEukPath[i]] for i in list_aa]
    parse.parrallelize(parse.dmdParse,to_leuk)
    print(f'step 9/9 :  ok parse dmd \n\n', file=open(log_f, 'a'))
    t1_stop = process_time()

#    t1_start = process_time()
 #   to_leuk = [[Models_df.aaPath[i], Models_df.dmdPath[i], Models_df.aaEukPath[i]] for i in list_aa]
  #  to_leuk=[to_leuk,args.lineage]
   # parse.dmdParse(to_leuk)
    #t1_stop = process_time()

    print('END !!!', file=open(log_f, 'a'))
    #     # # extract prodigal prot.
    #         # print(parse.extProd,Models_df['repo['aaEukPath]'])
    #         # print(parse.extProd,Models_df['repo['aaEukPath]'].values.tolist())
#     #         # parse.parrallelize(parse.extProd,Models_df['repo['aaEukPath]'].values.tolist())
#     #         #
#     #
#
#     T_stop = process_time()
#     print(f'\n\n---\nElapsed time: {t1_start} -> {t1_stop}', file=open(log_f, 'a'))
#     print(f'Elapsed time during the whole program in seconds:{T_stop}-{T_start}', file=open(log_f, 'a'))
#

### TEST
# print(gestError.krakenModels_dico)
