#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8   
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=20G
#$ -N clonalplots_VENCHI

conda activate r_4

echo "Started running clone plots for ${SGE_TASK_ID}"
Rscript /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code/Clonal_plots_paper.r $SGE_TASK_ID
echo "Completed running clones plots"
