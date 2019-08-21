import numpy as np
from pileup import load_allele_counts
from Bio import SeqIO
from Bio.Alphabet import generic_dna, Gapped
from Bio.Seq import Seq
import os
import glob


# Parse directories for which experiment and which vial
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Choose morbidostat experiment',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--exp_id',
                        help='ID of morbidostat experiment')
    parser.add_argument('--vial_nr',
                        help='Vial number of interest e.g. vial_1')
    args = parser.parse_args()
    exp_id = args.exp_id
    vial = args.vial_nr


base_dir = '/scicore/home/neher/ulreri00/morb_exp/'+exp_id+'/'+vial+'/'
samples = sorted(glob.glob(base_dir+'sample_*'))
alphabet_str = 'ACGT-N'

# Store pile_ups in matrices
def build_meta():
    meta = dict()
    fasta_file = list(SeqIO.parse(base_dir+'sample_0/unicycler/assembly.fasta','fasta'))
    ref = sorted([np.array(seq) for seq in fasta_file], key=lambda x:len(x),reverse = True)
    meta['ref_seq'] = ref
    contig_sizes = [len(seq) for seq in ref]
    contigs = len(contig_sizes)
    pile_ups = []
    for c in contig_sizes:
        pile_ups.append(np.zeros((len(samples),c,len(alphabet_str))))
    for sample_counter,sample in enumerate(samples):
        p = sample+'/pile_ups'
        pile_up = load_allele_counts(p)[0]
        pile_up.sort(key=lambda x:x[1].shape[-1],reverse=True)
        for ci, (cname,ac) in enumerate(pile_up):
            contig = len(pile_up[ci][1][-1])
            for pos in range(contig):
                pile_ups[ci][sample_counter][pos][:] = ac[:,pos]
    meta['pile_ups'] = pile_ups

    # Calculate attributes of pile_ups
    coverage = [np.repeat(pile_ups[c].sum(axis=2)[:,:,np.newaxis],6,axis=2) for c in range(contigs)]
    meta['coverage'] = coverage
    base_frequency = [pile_ups[c]/coverage[c] for c in range(contigs)]
    meta['base_freq'] = base_frequency
    base_index = [pile_ups[c].argmax(axis=2) for c in range(contigs)]
    meta['base_index'] = base_index

    # Find SNPs with a certain coverage and frequency
    filtred_pos = [[] for x in xrange(len(base_index))] 
    for contig_counter,contig in enumerate(base_index):
        for pos in range(len(contig[0])):
            bases = list(contig[:,pos])
            if min(bases) != max(bases) and min([c[0] for c in coverage[contig_counter][:,pos]])>30 and max([max(f) for f in base_frequency[contig_counter][:,pos]])>0.8:#np.amax(base_frequency[contig_counter][:,pos])>0.8:
                diff = min(bases, key=lambda x: (bases.count(x), bases[::-1].index(x)))
                diff_samples = [i for i,x in enumerate(bases) if x == diff]  
                filtred_pos[contig_counter].append(pos)
    meta['filtred_pos'] = filtred_pos
    return meta

def get_altered_genes(meta):
    gb_file = base_dir+'sample_0/annotations/prokka.gbk'
    genes = dict()
    not_annotated = dict()
    prev_pos = 0
    for contig_counter,gb_record in enumerate(SeqIO.parse(gb_file,'genbank')):
        for pos in meta['filtred_pos'][contig_counter]: 
            attributes = dict()       
            for feature in gb_record.features:
                if pos in feature and feature.type == 'CDS':
                    try:
                        genes[feature.qualifiers['product'][0]]['snp_pos'][2].append(pos)
                    except KeyError:
                        attributes['snp_pos'] = [contig_counter,feature.location.strand,[pos],[feature.location.start.real,feature.location.end.real]] 
                        attributes['nc_seq'],attributes['aa_seq'] = get_seq(meta['base_index'],attributes['snp_pos'])
                        try:
                            attributes['gene'] = feature.qualifiers['gene'][0]
                        except KeyError:
                            pass
                        genes[feature.qualifiers['product'][0]] = attributes
                else:
                    pass
    return genes 

def get_not_annotated(meta,genes):
    annotated = []
    all_snps = []
    for key in genes.keys():
        for pos in genes[key]['snp_pos'][2]:
            annotated.append((genes[key]['snp_pos'][0],pos))
    for contig_counter,contig in enumerate(meta['filtred_pos']):
        for pos in contig:
            all_snps.append((contig_counter,pos))
    not_ann = list(set(all_snps).difference(set(annotated)))

    not_annotated = dict()
    for (c,pos) in not_ann:
        not_annotated[pos] = [c,pos]
    return not_annotated


def get_seq(base_index,position):
    strand = position[1]
    contig = position[0]
    gene_pos = position[-1]
    nc_seq = dict()
    aa_seq = dict()
    for s, sample in enumerate(samples):
        if strand ==1:
            seq = Seq(''.join([alphabet_str[b] for b in base_index[contig][s][gene_pos[0]:gene_pos[1]]]))
            nc_seq[sample.split('/')[-1]] = seq
            aa_seq[sample.split('/')[-1]] = Seq(seq.tostring().replace('-','')).translate()
        if strand ==-1:
            seq = Seq(''.join([alphabet_str[b] for b in base_index[contig][s][gene_pos[0]:gene_pos[1]]])).reverse_complement()
            nc_seq[sample.split('/')[-1]] = seq
            aa_seq[sample.split('/')[-1]] = Seq(seq.tostring().replace('-','')).translate()
    return nc_seq,aa_seq

def generate_snp_csv(meta,genes):
    not_annotated = get_not_annotated(meta,genes)
    snp_counter = 0
    for key in genes:
        snp_counter += 1
        for pos in genes[key]['snp_pos'][2]:
            contig = genes[key]['snp_pos'][0]
            row = [str(snp_counter),str(contig),str(pos)]
            bases = [alphabet_str[b] for b in meta['base_index'][contig][:,pos]]
            for b in bases:
                row.append(b)
            print('\t'.join(row))

    for not_ann in not_annotated.values():
            snp_counter += 1
            row = [str(snp_counter),str(not_ann[0]),str(not_ann[1])]
            contig = not_ann[0]
            pos = not_ann[1]
            bases = [alphabet_str[b] for b in meta['base_index'][contig][:,pos]]
            for b in bases:
                row.append(b)
            print('\t'.join(row))

def generate_annotation_csv(meta,genes):
    for snp_counter,key in enumerate(genes):
        snp_counter+=1
        sample_len = len(genes[key]['aa_seq'])
        aa_len = len(genes[key]['aa_seq']['sample_0'])
        aa_matrix = np.chararray((sample_len,aa_len))
        for sample_index,sample in sorted(enumerate(genes[key]['aa_seq']),reverse=True):
            for base_index,base in enumerate(genes[key]['aa_seq'][sample]):
                aa_matrix[sample_index][base_index] = base
        for matrix_index in range(len(aa_matrix[0])):
            if len(set(aa_matrix[:,matrix_index])) > 1:
                aa = list(aa_matrix[:,matrix_index])
                aa = '\t'.join(aa)     
                try:
                    gene = genes[key]['gene']
                except KeyError:
                    gene = ''
                row = [str(snp_counter),gene,key,str(matrix_index),aa]
                print('\t'.join(row))      
