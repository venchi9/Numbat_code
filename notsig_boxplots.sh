#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8   
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=20G
#$ -N notisg_boxplots_VENCHI

conda activate R_4_2

echo "Started running notsig boxplots for ${SGE_TASK_ID}"
Rscript /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code/create_boxplots_notsig.r $SGE_TASK_ID
echo "Completed running notsig boxplots"
