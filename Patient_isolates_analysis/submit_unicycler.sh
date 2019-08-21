#!/bin/sh
#
# Reserve 16 CPUs for this job
#
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=1-0:00:00
#SBATCH --qos=1week
#



#move to current working directory
ml Unicycler/0.4.6-foss-2016b-Python-3.5.2
ml Pilon/1.22-Java-1.8.0_92
java -Xmx32G -jar $EBROOTPILON/pilon.jar
cd $SLURM_SUBMIT_DIR
unicycler -t 16 --short1 $1 --short2 $2 --long $3 --out $4 --spades_path /scicore/home/neher/GROUP/bin/anaconda2/bin/spades.py --racon_path /scicore/home/neher/GROUP/bin/anaconda2/bin/racon --pilon_path $EBROOTPILON/pilon.jar


