import glob
import os
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from shutil import copyfile
import pysam


base_dir = '/scicore/home/neher/ulreri00/morb_exp'+'/'
samples = sorted(glob.glob(base_dir+'0*/vial_*/sample*'))


def classify_reads():
    with open(base_dir+'read_counts.csv','a') as outfile:
        for sample in samples: 
            bacillus_bam = sample+'/mapped_to_bacillus.bam'
            ecoli_bam = sample+'/mapped_to_ecoli.bam'
            total_reads = pysam.view('-c',bacillus_bam).split('\n')[0]
            mapped_e_coli_rads = pysam.view('-c','-F','4',ecoli_bam).split('\n')[0]
            mapped_bacillus_reads = pysam.view('-c','-F','4',bacillus_bam).split('\n')[0]
            print(sample+'\t'+str(total_reads)+'\t'+mapped_bacillus_reads+'\t'+mapped_e_coli_rads+'\n')
            outfile.write(sample+'\t'+str(total_reads)+'\t'+mapped_bacillus_reads+'\t'+mapped_e_coli_rads+'\n')

                    