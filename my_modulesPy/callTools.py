#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import pathlib
from modules import parse

#run augustus
def augustus(lofl): ## ajouter singularity avec variable sing et sing2
    model,fastaN, gffPath=lofl
    if pathlib.Path(fastaN+'.fna').exists ():
        str1 = f'augustus --uniqueGeneId=true --gff3=on --species={model} {fastaN}.fna  > {gffPath}.gff'
        content = os.popen(str1).read()


# dev!! a verifier
def augustusLoop2(dt):   #fastaRepo,dico_contigTaxon,fastaOutRepo
    job_list = []
    for aug_model in dt.index.to_list():
        t_list = [aug_model,dt.fnaPath[aug_model],dt.gffPath[aug_model]]
        job_list.append(t_list)
    parse.parrallelize(augustus, job_list)


def addH(f,i): # add metagid in fasta header
    str1 = "./addH.sh "+f+' '+i
    content = os.popen(str1).read()

def diamond(l): ## ajouter singularity avec variable sing et sing2
    aa=f'{l[0]}.aa'
    dmd=f'{l[1]}.dmd'
    db='/lerins/hub/DB/NR/NR_diamond/NR_2020_01_diamond'
    if not pathlib.Path(db).exists():
        str1 = f'diamond blastp -p 60   -e 0.000001 --outfmt 102 -d {db} -q {aa} -o {dmd}'
        content = os.popen(str1).read()
    else:
        print(f'error : {db=} doesn\'t exist ')
