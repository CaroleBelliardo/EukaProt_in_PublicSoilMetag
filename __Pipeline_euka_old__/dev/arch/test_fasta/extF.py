#!/usr/bin/env python3
# Execute by: python seq_extractor.py readsList fastafile outputfile

from Bio import SeqIO
import sys

readsList = open(sys.argv[1], 'r')
fastafile = sys.argv[2]
outputfile = open(sys.argv[3], 'w')

wanted = set()
with readsList as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fastafile),'fasta')

with outputfile as i:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], i, "fasta")

readsList.close()
outputfile.close()
