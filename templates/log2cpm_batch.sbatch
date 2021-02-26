#!/bin/bash

#SBATCH --array=1-56%30
#SBATCH --mem=3G
#SBATCH -J log2cpm
#SBATCH -o /scratch/mblab/$USER/rnaseq_pipeline/sbatch_log/log2cpm.out

ml rnaseq_pipeline

read count_file < <(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/mblab/chasem/rnaseq_pipeline/job_scripts/induction_set_countfile_lookup.txt )

organism=KN99

log2_cpm.R -r ${count_file} -g ${organism} 

