#!/bin/sh
#
# Reserve 16 CPUs for this job
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=05:29:00
#
source activate panX
#move to current working directory
cd $SLURM_SUBMIT_DIR
./panX.py -fn data/e_coli -sl e_coli -t 16 -iba -mi data/e_coli/metainfo_esbl.tsv -mtf data/e_coli/metadata_config_esbl.tsv -st 9 10 11 > e_coli.log 2> e_coli.err 