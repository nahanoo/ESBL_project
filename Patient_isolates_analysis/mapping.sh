#!/bin/sh
#
# Reserve 16 CPUs for this job
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=00:29:00
#
#move to current working directory
cd $SLURM_SUBMIT_DIR
ml BWA
ml SAMtools/1.7-goolf-1.7.20
bwa mem $1 $2 $3 | samtools view - -Sb | samtools sort - -o $4