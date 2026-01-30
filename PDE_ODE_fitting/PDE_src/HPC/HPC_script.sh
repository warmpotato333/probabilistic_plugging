#!/bin/bash
#SBATCH --job-name kumo5_la12_Bo1_REDO_job
#SBATCH --output Matlab_kumo5_la12_Bo1_job.log
#SBATCH --error Matlab_kumo5_la12_Bo1_job.log
#SBATCH --cpus-per-task=128
#SBATCH --mem 256G
 
module load matlab/R2024a
matlab -nodisplay -nojvm  < kumo6_LW_Flow_SweepA_la12_Bo1.m