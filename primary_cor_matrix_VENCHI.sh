#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 8   
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=20G
#$ -N primary_cor_matrix_VENCHI

conda activate R_4_2

echo "Started running primary cor_matrix for ${SGE_TASK_ID}"
Rscript /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code/Corr_matrix_mean_primary.r $SGE_TASK_ID
echo "Completed running primary cor_matrix Objects"
