# EukaProt_in_PublicSoilMetag

Bioinformatic python3 scripts to predict eukaryotic proteins in Metagenome-Assembled Genomes (MAGs) data. 


Dependances: 
  Python3 package : 
 numpy, 
   - pandas,
   - biopython
   - pyfasta
   - time 
   - argparse
   - pprint
    

  This pipeline use several tool : 
   - KRAKEN2 (contig taxonomic classification), 
   - AUGUSTUS (protein prediction), 
   - DIAMOND (protein taxonomic annotation),
   - BUSCO (Protein prediction quality). 
