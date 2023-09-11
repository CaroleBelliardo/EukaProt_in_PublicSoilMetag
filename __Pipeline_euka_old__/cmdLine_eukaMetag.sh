# --Contigs fasta files download
# use API
# --run Kraken
## DB format
# singularity kraken require blast for dust...and seq..
#------------------------------------------------------------------------------------------------------------
#!/bin/bash

SING_IMG=/bighub/hub/people/carole/work-bighub/sing-image/MetagTools_kraken_dmd_blast.sif
SING2="singularity exec --bind /bighub/hub:/bighub/hub:rw --bind /work/cbelliardo:/work/cbelliardo:rw"
#SING2="singularity exec"

DBNAME=/bighub/hub/DB/refseq_kraken

$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-taxonomy --db $DBNAME --threads 70

$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library bacteria --db $DBNAME --threads 70
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library archaea --db $DBNAME --threads 70
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library plasmid --db $DBNAME --threads 70 
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library viral --db $DBNAME --threads 70
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library human --db $DBNAME --threads 70
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library fungi --db $DBNAME --threads 70
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library plant --db $DBNAME --threads 70
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library protozoa --db $DBNAME --threads 70
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --download-library nt --db $DBNAME --threads 70

$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-build --build --db $DBNAME --threads 70


# Faire un krona de la db pour voir si compo = biaise results
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2-inspect --db $DBNAME > k


# kraken to krona
out=All_0.1.krak
cat $1 > $out
cut -f 3 $out | sort | uniq -c | sed 's/^ *//'  | sed  's/ /\t/'  > ${out}_cut
awk '{ FS = OFS = "\t" } NR==FNR {h[$1] = $2 ; next} {print $0,h[$2]}' $db_lineage ${out}_cut | cut -f 1,3 > ${out}_cut_lineage
sed -i 's/\t;/\t/' ${out}_cut_lineage
sed -i 's/;$//' ${out}_cut_lineage
sed -i 's/;/\t/g' ${out}_cut_lineage
grep -vi 'root' ${out}_cut_lineage > ${out}_cut_lineage_class 
ktImportText ${out}_cut_lineage_class -o ${out}_cut_lineage_class.html



# run kraken
$SING2 $SING_IMG /home/tools/kraken2-master/scripts/kraken2 --confidence 0.1 --db $DBNAME --threads 70 AllrenameHeader_data_*.fna


#-- RunPipline
SING_IMG='/bighub/hub/people/carole/work-bighub/sing-image/MetagTools_kraken_dmd.sif'
SING2="singularity exec --bind /bighub/hub:/bighub/hub:rw --bind /work/cbelliardo:/work/cbelliardo:rw"

wd='/bighub/hub/people/carole/work-bighub/Scripts_git/Pipeline_euka/dev'
fnaRepo='/work/cbelliardo/1-seuil_cut-Metag/Faa_1000l-3g_paquets'
modelAugustus='/bighub/hub/people/carole/work-bighub/Scripts_git/Pipeline_euka/dev/tx_sp_mod_nods.txt'
kraken_euka='/work/cbelliardo/3-PipelineProtEuka/V2/1-kraken/test_db_confidence/AllrenameHeader_data_0.1.krak'
tool='/bighub/hub/people/carole/work-bighub/Scripts_git/Pipeline_euka/dev/Branch_pipEuka.py'
output=$1

$SING2 $SING_IMG python3 -u $tool -f $fnaRepo -mt $modelAugustus -k $kraken_euka -o $output >> ${1}.log



