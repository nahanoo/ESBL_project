#!/bin/sh
#
# Reserve 1 CPUs for this job
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=00:59:00
#


#move to current working directory
cd $SLURM_SUBMIT_DIR

ml Trim_Galore
echo "trim_galore $@"
trim_galore $@

