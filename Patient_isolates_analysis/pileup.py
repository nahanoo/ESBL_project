#!/usr/bin/env python2
import numpy as np
from collections import defaultdict


def load_allele_counts(dirname, suffix=''):
    import cPickle, gzip, glob
    dirname = dirname.rstrip('/')+'/'
    tmp_ac = {}
    ac_flist = glob.glob(dirname+'*_allele_counts' + suffix + '.npz')
    for fname in ac_flist:
        #print("reading",fname)
        tmp = '_allele_counts' + suffix + '.npz'
        refname = fname.split('/')[-1][:-len(tmp)]
        print(fname)
        tmp_ac[refname] = np.load(fname).items()[0][1]

    ins_flist = glob.glob(dirname+'*insertions' + suffix + '.pkl.gz')
    tmp_ins = {}
    for fname in ins_flist:
        print("reading",fname)
        tmp = '_insertions' + suffix + '.pkl.gz'
        refname = fname.split('/')[-1][:-len(tmp)]
        with gzip.open(fname) as fh:
            tmp_ins[refname] = cPickle.load(fh)

    ac = []
    for refname in tmp_ac:
        ac.append((refname, tmp_ac[refname]))

    # ins = []
    # for refname in tmp_ins:
    #     ins.append((refname, tmp_ins[refname]))

    return ac, tmp_ins


def dump_allele_counts(dirname, ac, suffix=''):
    import cPickle, gzip, os
    dirname = dirname.rstrip('/')+'/'
    if not os.path.isdir(dirname):
        print("creating directory", dirname)
        try:
            os.mkdir(dirname)
        except:
            print("creating directory failed", dirname)

    for refname, ac_array, insertions in ac:
        print(refname)
        outname = refname.split("|")[-1]
        np.savez_compressed(dirname + outname+'_allele_counts' + suffix + '.npz', ac_array)
        #with gzip.open(dirname + outname+'_insertions' + suffix + '.pkl.gz','w') as outfile:
        #    cPickle.dump({k:dict(v) for k,v in insertions.iteritems()}, outfile)


def returnPileupStructure( N, QC, pairedEnd ):
    if (not QC):
        return np.zeros((6,N))
    else:
        if (pairedEnd):
            return np.zeros((2,2,6,N))
        else:
            return np.zeros((2,6,N))

def returnInsertionStructure( QC, paired ):
    if QC:
        if paired:
            return defaultdict(lambda: defaultdict(lambda: np.zeros((2,2),int)))
        else:
            return defaultdict(lambda: defaultdict(lambda: np.zeros(2,int)))
    else:
        return defaultdict(lambda: defaultdict(lambda: np.zeros(1,int)))

def performPileUp( bamFile, QC=False, pairedEnd=False, qual_min=0, max_reads=-1):
    import pysam
    alphabet = np.fromstring('ACGT-N', dtype='S1')
    align = pysam.AlignmentFile(bamFile,"rb")
    # Allocate space for pileup
    totCov = []
    for nContig in xrange(align.nreferences):
        totCov.append( [align.getrname(nContig),
                        returnPileupStructure(align.lengths[nContig],QC,pairedEnd),
                        returnInsertionStructure(QC,pairedEnd)] )

    # Iterate over reads in BAM file
    nUnmapped = 0
    nReads=0
    for n, read in enumerate(align):
        nReads+=1
        if (max_reads>0 and nReads>max_reads):
            break
        if read.seq is None:
            continue
        if read.is_unmapped:
            nUnmapped+=1
            continue

        ## shorthand for the relevant data structures
        if (QC):
            rev = int(read.is_reverse)
            if (pairedEnd):
                r1 = int(read.is_read2)
                counts = totCov[read.rname][1][r1,rev]
            else:
                counts = totCov[read.rname][1][rev]
        else:
            counts = totCov[read.rname][1]

        insertion = totCov[read.rname][2]

        # convert sequence and quality to arrays
        seq = np.fromstring(read.seq,'S1')
        if (read.qual is not None):
            qual = np.fromstring(read.qual,np.int8) - 33
        else:
            qual = 60*np.ones(len(seq))
        pos = read.pos

        # Iterate over block types within the CIGAR code
        for ic, (block_type, block_len) in enumerate(read.cigar):
            if (block_type == 4): # Soft clip
                seq = seq[block_len:]
                qual = qual[block_len:]

            if (block_type == 0): # Normal
                seqb = seq[:block_len]
                qualb = qual[:block_len]
                for j, a in enumerate(alphabet): # Increment nucleotide counts from read
                    posa = ((seqb==a) & (qualb >= qual_min)).nonzero()[0]
                    if len(posa):
                        counts[j,pos+posa] += 1

                if ic != len(read.cigar) - 1: # Move along sequence
                    seq = seq[block_len:]
                    qual = qual[block_len:]
                    pos += block_len

            elif (block_type == 2): # Deletion
                counts[4,pos:pos+block_len] += 1
                pos += block_len #Advance reference but not read

            elif (block_type == 1): # Insertion
                seqb = seq[:block_len]
                qualb = qual[:block_len]
                # Accept only high-quality inserts
                if (qualb >= qual_min).all():
                    if (QC):
                        if pairedEnd:
                            insertion[pos][seqb.tostring()][r1, rev] += 1
                        else:
                            insertion[pos][seqb.tostring()][rev] += 1
                    else:
                        insertion[pos][seqb.tostring()] += 1

                # Chop off seq, but not pos
                if ic != len(read.cigar) - 1:
                    seq = seq[block_len:]
                    qual = qual[block_len:]

    fracUnmapped = (nUnmapped) / (1. * n)

    return totCov, fracUnmapped

if __name__ == '__main__':
    from Bio import SeqIO
    import argparse
    parser = argparse.ArgumentParser(description='create allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bam-file',
                        help='bam file to pile up')

    parser.add_argument('--out-dir',
                        help='directory to save results')
    args = parser.parse_args()

    #fname = '../ulreri00/ESBL_project/Patients/Patient02/Sample1/samples_polished.bam'
    fname = args.bam_file

    pileup, f = performPileUp(fname)


    dump_allele_counts(args.out_dir,pileup)