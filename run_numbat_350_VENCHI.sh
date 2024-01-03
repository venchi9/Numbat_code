#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 30   
#$ -cwd
#$ -V
#$ -q short.q
#$ -r yes
#$ -l mem_requested=20G
#$ -N numbat_VENCHI

conda activate r_4

echo "Started running numbat for ${SGE_TASK_ID}"
Rscript /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/code/run-numbat_350_VC.R $SGE_TASK_ID
echo "Completed running numbat"
