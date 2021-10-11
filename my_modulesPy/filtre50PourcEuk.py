#!/usr/bin/python3
# -*- coding: utf-8 -*-

##--- Headercd 
import pprint   # debugg -- affichage Dump des structures de données
import os       # manip systeme
import sys      # gestion arguments et opt

## init
Fi=sys.argv[1]
outF=sys.argv[1]+'min50PEuk'

File = open(Fi, "r")   ## travail sur version triée de res. de blast
outFile = open(outF, "w")  

dico = dict()
for ligne in File: # -- Parcourt du fichier d'entrée ( blast trié )
	ligne=ligne.rstrip()    # supprime les sauts de lignes
	l=ligne.split("\t")     # coupe la ligne en liste au niveau des tab
	if (float(l[2]) >= 0.5*int(l[0])) :
		outFile.write(l[1]+"\n")


File.close()
outFile.close()
