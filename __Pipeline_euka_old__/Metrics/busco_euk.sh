#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --cpus-per-task=4     # number of CPU per task
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=32G   # memory per Nodes
#SBATCH -J "busco"   # job name
#SBATCH --mail-user=carole.belliardo@inra.fr   # email address
#SBATCH --mail-type=ALL
#SBATCH --gid=bioinfo
#SBATCH -e /work/cbelliardo/zslurm-jobs/slurm-%j.err
#SBATCH -o /work/cbelliardo/zslurm-jobs/slurm-%j.out
#SBATCH -p all

SING_IMG=/bighub/hub/people/carole/work-bighub/sing-image/busco_v4.0.5_cv1.sif
SING2="singularity exec --bind /bighub/hub:/bighub/hub:rw --bind /work/cbelliardo:/work/cbelliardo:rw"

cd $1
# grab our filename from a directory listing
FILES=($(ls -1 *.aa))
MODELS=($(ls -1 *.aa | cut -d"_" -f1))
#METAGS=($(ls -1 *.aa | cut -d"." -f2))
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}
#MODEL=${MODELS[$SLURM_ARRAY_TASK_ID]}


#busco --config /bighub/hub/people/carole/work-bighub/sing-image/busco_downloads
lc=/bighub/hub/people/carole/work-bighub/sing-image/busco_downloads/lineages
$SING2 $SING_IMG busco -i $FILENAME -c 4 -o ${FILENAME}.busco -l ${lc}/eukaryota_odb10 -m protein #--augustus_species $MODEL 
# $FILENAME > ${model}_${PREF}
#alveolata_odb10






#lancement du script:
#while read line ; do ln -s ../$line.fa . ; done < LS-split-1689
#ln -s ../*fa .

# get count of files in this directory
#NUMFILES=$(ls -1 *.fa | wc -l)
# subtract 1 as we have to use zero-based indexing (first element is 0)
#ZBNUMFILES=$(($NUMFILES - 1))


#sbatch -p all --array=0-$ZBNUMFILES%120 augustus.sh

