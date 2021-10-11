#!/bin/bash

f=$1
db_lineage='/lerins/hub/DB/TAXONOMY/online/taxonomy_lineage.txt'
awk '{ FS = OFS = "\t" } NR==FNR {h[$1] = $2 ; next} {print $0,h[$2]}' $db_lineage $f > $f.lineage

