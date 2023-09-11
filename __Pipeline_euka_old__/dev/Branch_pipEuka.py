#!/usr/bin/python3
# -*- coding: utf-8 -*-

import pprint  # debugg -- affichage Dump des structures de données
import os  # manip systeme
import pathlib  # test file
import sys  # gestion arguments et opt
import argparse
import numpy as np
import pandas as pd
import re
from modules import gestError
from modules import parse
from modules import callTools
from time import process_time


# Options
parser = argparse.ArgumentParser(
    description='Need PYTHON 3 version with packages numpy, pandas, biopython install with :pip install biopython; '
                'pip install numpy; pip install pandas; Ex. command line :  python3 Pipeline_euk_v2.py -fna '
                'fasta_test -mt tx_sp_mod_nods.txt -k krak.test -t Mucoromycota -o test ')
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
# # Main
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    with open(log_f, 'w') as log:
        # -- Kraken to fna datasets
        t1_start = process_time()  # need process_time
        gestError.file_exist(args.modelTable)  # -- import augusutus_modelsTable
        Models_df = pd.read_table(args.modelTable, sep='\t')
        Models_df = Models_df.set_index('model')
        t1_stop = process_time()
        log.write(f'step1 : ok index  {args.modelTable}\n      time krakenLineage: {t1_stop - t1_start}')  # rownames model table

        # -- List euka contigs
        t1_start = process_time()
        krakenModels_dico, message = gestError.krakenModels(args.kraken, args.lineage, tmp_r, Models_df, args.eukTaxid)
        t1_stop = process_time()
        log.write(f'step2 : {message} ok kraken dico  {args.kraken}_lineage\n'); log.write(f'     time krakenLineage: {t1_stop - t1_start}  \n'); log.write(f'    {len(args.kraken)=} \n')
        print(krakenModels_dico)

        # -- parse kraken output



    # => {model: ['Ga0272421_1202915_3300031909'], 10090: ['Ga0272421_1213909_3300031909'] ...

    # #-- Paths
    #     fastaOutRepo=args.output+'/FNA_krakenEuka-contigs/'; gestError.mkdir_exist(fastaOutRepo) ;print('step8 : ok CreateRepo kraken_fna')
    #     # creation sous dossier
    #     dtdt= Models_df.deeper_taxon.str.replace(' ','_')
    #     Models_df['fnaPath'] = fastaOutRepo+dtdt.astype(str); Models_df.fnaPath.apply(gestError.mkdir_exist);print('step9 : ok CreateSousRepo kraken_fna')
    #     GffOutRepo=args.output+'/GFF_krakenEuka-contigs/'; gestError.mkdir_exist(GffOutRepo) ;print('step7 : ok CreateRepo kraken_gff')
    #     Models_df['gffPath'] = GffOutRepo + dtdt.astype(str) +'.gff'
    #     BEDOutRepo=args.output+'/BED_krakenEuka-contigs/'; gestError.mkdir_exist(BEDOutRepo) ;print('step9: ok CreateRepo BED')
    #     Models_df['bedPath'] = BEDOutRepo + dtdt.astype(str) +'.bed'
    #     AAOutRepo=args.output+'/AA_krakenEuka-contigs/'; gestError.mkdir_exist(AAOutRepo) ;print('step10 : ok CreateRepo AA')
    #     Models_df['aaPath'] = AAOutRepo + dtdt.astype(str) +'.aa'
    #     DMDOutRepo=args.output+'/DMD_krakenEuka-contigs/'; gestError.mkdir_exist(DMDOutRepo) ;print('step12: ok CreateRepo DMD')
    #     Models_df['dmdPath'] = DMDOutRepo + dtdt.astype(str) +'.dmd'
    #     AaEukPath=args.output+'/AAeuk_krakenEuka-contigs/'; gestError.mkdir_exist(AaEukPath) ;print('step10 : ok CreateRepo AA')
    #     Models_df['aaEukPath'] = AaEukPath + dtdt.astype(str) +'.aa'
    # # #-- extract seq euka
    #     # parse.extractSeqRun(args.fna,Dico_contigTaxon,fastaOutRepo,Models_df[['fnaPath','deeper_taxon']]) ;print('step10 : ok extract euka seq ')
    # # #-- concatenate fasta: one file by model
    #     parse.parrallelize(parse.concatFile,Models_df['fnaPath'].tolist()); print('step11 : ok concat fasta files')
    # # # -- dev ----***
    # # #-- check file size
    # #     ## MODIF. => Pas de creation repo mais split file dans Repo courantet ajout a la liste des fichier;
    # #     # Pour chaque modekl, PAth variable == List!!
    # #     # ajout iteration sur liste repo à toutes les etapes !!
    # #         # parse.parrallelize(parse.checkFile_split,Models_df['fnaPath'].tolist()); print('step7.5 : ok check fasta files size')
    # # # -- dev end ----***
    # #-- run augustus
    #     callTools.augustusLoop2(Models_df)
    # # -- parse augfiles
    #     to_lgff=pd.DataFrame(Models_df, columns= ['gffPath','bedPath', 'aaPath']).values.tolist()
    #     to_lgff=[GffOutRepo,BEDOutRepo,AAOutRepo]
    #     parse.gffParse2(to_lgff);print('step11 : ok parse Gff')
    #  #-- run diamond
    #     to_ldmd=pd.DataFrame(Models_df, columns= ['aaPath', 'dmdPath']).values.tolist()
    #     for ilist in to_ldmd : callTools.diamond(ilist)
    #     to_leuk=pd.DataFrame(Models_df, columns= ['aaPath', 'dmdPath','aaEukPath']).values.tolist()
    #     parse.parrallelize(parse.dmdParse,to_leuk)
    # # # extract prodigal prot.
    # #         # print(parse.extProd,Models_df['aaEukPath'])
    # #         # print(parse.extProd,Models_df['aaEukPath'].values.tolist())
    # #         # parse.parrallelize(parse.extProd,Models_df['aaEukPath'].values.tolist())
    # #         #
    # #
    # t1_stop = process_time()
    # print("\n\n---\nElapsed time:", t1_start, '->',t1_stop)
    # print("Elapsed time during the whole program in seconds:",
    #                                          t1_stop-t1_start)
