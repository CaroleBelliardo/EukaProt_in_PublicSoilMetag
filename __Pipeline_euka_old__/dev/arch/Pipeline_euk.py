#!/usr/bin/python3
# -*- coding: utf-8 -*-

##--- Header
import pprint   # debugg -- affichage Dump des structures de données
import os       # manip systeme
import sys      # gestion arguments et opt
import argparse
import numpy as np
import pandas as pd

# options
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-mt', '--modelTable', type=str, help='modele - species (.tab file) \n')
parser.add_argument('-k', '--kraken', type=str, help='Kraken output (.tab file) \n')
parser.add_argument('-t', '--taxon', type=str, help='taxon selected for only one taxon (.tab file) \n')
parser.add_argument('-ps', '--printS', type=bool, nargs='?', default=False, help='modele - species (.tab file) \n')
args = parser.parse_args()


## init model Augustus table
dt=pd.read_table(args.modelTable, sep='\t') # import model table 
dt = dt.set_index('species') # rownames

def indique_mon_age():
	return 30


#if args.printS != False: #print liste espèces modèles
	





#Filtre contig Kraken selon esp
#Open kraken file
#os.system ("bash -c 'sed '1s/^/stat\tcontig\tf\ttaxid\tlist\tlineage\n/' krak.test > krak.test_h'")

str1 = "sed '1s/^/stat\\tcontig\\tf\\ttaxid\\tlist\\tlineage\\n/' "+args.kraken+">"+args.kraken+"_h"
content = os.popen(str1).read()
dk = pd.read_table(args.kraken+"_h", skipinitialspace=True, usecols=['contig', 'lineage'], sep='\t')











#outFile = open(outF, "w")  
#MT = open(args.modelTable, "r")   ## travail sur version triée de res. de blast
#dico = dict()
#for ligne in MT: # -- Parcourt du fichier d'entrée ( blast trié )
#	ligne=ligne.rstrip()    # sudtrime les sauts de lignes
#	l=ligne.split("\t")     # coupe la ligne en liste au niveau des tab
#	print(l)
#File.close()
#outFile.close()

## Main 
# if species exit alors ex pour une esp
# sinon all speciesmodelTable

if __name__ == "__main__":
	if args.printS != False: #print liste espèces modèles
		pd.set_option('display.max_rows', dt.shape[0]+1)
		print(dt)
### 
	if args.taxon == None :
		print('No specific taxon, analyse of all taxon*')
		# dérouler boucle sur tous les taxons
	else :
		sst=dt.loc[dt['deeper_taxon'] == args.taxon]
		if sst.empty:
			print(args.taxon+' Not in taxon list, use -ps postion to see all available taxon_node !!!')		
		else:
			print(args.taxon+' in list of taxon***')    
