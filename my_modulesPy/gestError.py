#!/usr/bin/python3
# -*- coding: utf-8 -*-
import pathlib  # test file
import os
import sys
import pandas as pd

def testFun2(k,v_set):
    dic[k]=v_set

def mkdir_exist(repoP):  # creation repertoire
    if not os.path.exists(repoP):
        os.makedirs(repoP)


def file_exist(filepath):  # test if file existe
    file_p = pathlib.Path(filepath)
    if not file_p.exists():
        sys.exit("Error message:\n\n!!! " + filepath + "file not exist !!!\n\n")
    return file_p


def file_exit_func(filepath, func):  # test if file existe
    file_p = pathlib.Path(filepath)
    if file_p.exists():
        print('  !! ' + filepath + ' ok')
        r = func(filepath)
    else:
        sys.exit("Error message:    \n\n!!! " + filepath + "file not exist !!!\n\n")
    return r


def KrakenLineage_exist(kraken_f, lineage_f, tmp_repo, log_f):  # -- add lineage to kraken C
    krakenLineage_f = tmp_repo + '/' + os.path.basename(kraken_f) + "_lineage.krak"
    print(f'    {krakenLineage_f} importing step :', file=open(log_f, 'a'))
    if not pathlib.Path(krakenLineage_f).exists():
        print(f'        {krakenLineage_f} not exist; creating... ', file=open(log_f, 'a'))
        kraken_df = pd.read_table(kraken_f, sep='\t', names=['classe', 'seqid', 'taxid', 'seqLen', 'LCAListe'])
        print(f'ok\n        {kraken_f}  filtering... ', file=open(log_f, 'a'))
        kraken_df = kraken_df[kraken_df.classe == 'C']
        print(f'ok\n        {kraken_f}  filtrering col...', file=open(log_f, 'a'))
        kraken_df[['seq', 'metag']] = kraken_df.seqid.str.rsplit('_', n=1, expand=True)
        print(f'ok\n        {lineage_f}  importing ... ', file=open(log_f, 'a'))
        lineage_df = pd.read_table(lineage_f, sep='\t', names=['taxid', 'lineage', 'nodes'])
        print(f'ok\n        {krakenLineage_f} and {lineage_f}   merging...  ', file=open(log_f, 'a'))
        kraken_lineage = pd.merge(kraken_df, lineage_df, how='left', on='taxid')
        del lineage_df
        print(f'ok\n        {krakenLineage_f}  printing... ', file=open(log_f, 'a'))
        kraken_lineage.to_csv(krakenLineage_f, index=False, sep='\t')
        print(f'        {krakenLineage_f}  printing ok !!! end step lineage paste kraken ', file=open(log_f, 'a'))
    else:
        print(f'        {krakenLineage_f} exist; importing... ', file=open(log_f, 'a'))
        kraken_lineage = pd.read_table(krakenLineage_f, sep='\t')
        print(f'        {krakenLineage_f}   imported ok ', file=open(log_f, 'a'))
    return kraken_lineage


# def format_kraken_dico(kraken_df, kraken_model_dico):
#     km_dico_formated = dict()
#     tmp = dict()
#     for model, contig_set in kraken_model_dico.items():
#         for m in set(kraken_df.metag[kraken_df.seqid.isin(contig_set)]):
#             contigs_m_set = set(kraken_df.seq[kraken_df.seqid.isin(contig_set)])
#             km_dico_formated = parse.appendDicoIndicoValue(km_dico_formated, m, model, contigs_m_set)
#     return km_dico_formated # ex:    km_dico={'fileB':{'rhizopus_oryzae': {'Ga0272421_100ustr', 'Ga0272421_100usta', 'Ga0272421_1000028', 'Ga0272421_100ustc', 'Ga0272421_1000035', 'Ga0272421_1000029'}}}
#
