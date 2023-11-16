#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8   
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=20G
#$ -N clone_heatmaps_VENCHI

conda activate r_4

echo "Started running heatmaps Objects for ${SGE_TASK_ID}"
Rscript /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code/Clone_heatmaps.r $SGE_TASK_ID
echo "Completed running heatmaps Objects"
