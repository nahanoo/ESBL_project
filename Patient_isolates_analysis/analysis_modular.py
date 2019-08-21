import numpy as np
from pileup import load_allele_counts
from Bio import SeqIO
import os
import glob

base_dir = os.getcwd()+'/Patients/'
patients = ['Patient33']
#,'Patient12','Patient25','Patient33'
alphabet_str = 'ACGT-N'
alphabet = np.fromstring(alphabet_str, dtype='S1')

def load_ref_genome(patient):
    # load fasta file from ref genome, patient number necessary E.g Patient01, also used for counting contigs
    cwd = base_dir+patient+'/'
    ref = list(SeqIO.parse(cwd+'Sample0/assembly.fasta','fasta'))
    ref_genome = sorted([np.array(seq) for seq in ref], key=lambda x:len(x))
    contig_sizes = []
    for rec in ref:
        contig_sizes.append(len(rec))
    contig_sizes = sorted(contig_sizes)
    return ref_genome,contig_sizes


def load_pile_ups(patient):
    # load pile ups from all samples, patient number necessary E.g Patient01
    # pile_ups are loaded into a list with all contigs, in each element of this list (contig), \
    # all samples are saved with the counts for the base. Alphabet is ACGT-N 
    contig_sizes = sorted(load_ref_genome(patient)[1],reverse = True)
    cwd = base_dir+patient+'/'

    sample_list = sorted(os.listdir(cwd))[1:] # exclude Ref_genome
    pile_ups = []
    for n in range(len(contig_sizes)):
        contigs = np.zeros((len(sample_list),contig_sizes[n],len(alphabet_str)))
        pile_ups.append(contigs)

    for sample in range(len(sample_list)):
        pile_up_name = cwd+sample_list[sample]+'/Pile_ups'
        pile_up = load_allele_counts(pile_up_name)[0]
        pile_up.sort(key=lambda x:x[1].shape[-1],reverse=True)

        for ci, (cname,ac) in enumerate(pile_up):
            length = len(pile_up[ci][1][-1])
            for pos in range(length):
                pile_ups[ci][sample][pos][:] = ac[:,pos]
                
    return pile_ups

def get_coverage(pile_ups):
    # returns coverage for each position in same data structure as the the pile ups
    coverage = list()
    for contig in range(len(pile_ups)):
        cov = pile_ups[contig].sum(axis=2)
        cov = np.repeat(cov[:,:,np.newaxis],6,axis=2)
        coverage.append(cov)
    return coverage

def get_base_frequency(pile_ups,coverage):
    # calculates base frequency
    base_freq = list()
    for contig in range(len(pile_ups)):
        freq = pile_ups[contig]/coverage[contig]
        base_freq.append(freq)
    return base_freq
        
def get_most_abundant_base(pile_ups):
    base_index = list()
    for contig in range(len(pile_ups)):
        base = pile_ups[contig].argmax(axis=2)
        #base = np.repeat(base[:,:,np.newaxis],6,axis=2)
        base_index.append(base)
    return base_index

def get_different_positions(base_index,pile_ups):
    # get different positions based on most abundant base
    diff_pos = [[] for x in xrange(len(pile_ups))]
    for contig in range(len(pile_ups)):
        for sample in range(len(pile_ups[0])):
            for sample_to_compare in range(len(pile_ups[0])):
                diff = np. where(base_index[contig][sample] != base_index[contig][sample_to_compare])[0]
                if len(diff)>0 and any(np.isin(diff,diff_pos[contig]))==False:
                    diff_pos[contig].append(diff)
    return diff_pos

def filter_diff_pos(pile_ups,coverage,base_freq,diff_pos,min_cov,min_base_freq):
    # filters different positions where coverage is below a certain coverage or base frequency
    filtred_pos = [[] for x in xrange(len(pile_ups))]
    for contig in range(len(pile_ups)):
        if diff_pos[contig]:
            for pos in diff_pos[contig][0]:
                if np.amax(coverage[contig][:,pos])>min_cov and np.amax(base_freq[contig][:,pos])>min_base_freq:
                    filtred_pos[contig].append(pos)
        return filtred_pos


def get_genes(patient,filtred_pos):
    # reads in the anntoated genbank file for the patient and goes through all the different positions
    # returns the gene name, the product and the EC_number 
    gb_file = base_dir+patient+'/Ref_genome/Annotations/prokka.gbk'
    # puts all the different positinos into one list
    genes = [] 
    prev_pos = 0
    del_counter = 0
    for contig_counter,gb_record in enumerate(SeqIO.parse(open(gb_file,'r'),'genbank')):
        for pos in filtred_pos[contig_counter]:                                                            
            for feature in gb_record.features:
                if pos in feature:
                    if prev_pos+1 !=  pos:
                        try:
                            location = [feature.location.start.real,feature.location.end.real]
                            gene_info = [pos,contig_counter,feature.qualifiers['gene'][0],feature.qualifiers['product'][0],feature.location.strand,location]
                            try:
                                if not genes[-1][2]==feature.qualifiers['gene'][0]:
                                    genes.append(gene_info)
                            except IndexError:
                                genes.append(gene_info)
                        except KeyError:
                            pass
                    else:
                        del_counter += 1              
            prev_pos = pos
    return genes

def translate_seq_from_pile_ups(patient,base_index,filtred_pos,gene_of_interest):
    from Bio import Seq
    
    # find the gene_of interest
    
    gene_location = gene_of_interest[-1]
    contig = gene_of_interest[1] 
    strand = gene_of_interest[-2]

    samples = sorted(os.listdir(base_dir+patient))[1:]
    nucleotide_sequences = {}
    keys = samples
    value = None
    

    # generate the nucleotide sequence dictionary and save the sequence
    for key in keys:
        nucleotide_sequences[key] = value

    for sample in samples:
        seq_trans = []
        sample_index = samples.index(sample)
        seq = base_index[contig][sample_index][gene_location[0]:gene_location[1]]
        for pos in seq:
            seq_trans.append(alphabet_str[pos])
        nucleotide_sequences[sample] = ''.join(seq_trans)

    
    # generate the amino acid sequence dictionary and save the sequence
    amino_acid_sequences = {}
    keys = samples
    value = None

    for key in keys:
        amino_acid_sequences[key] = value

    for sample in samples:
        seq = nucleotide_sequences[sample]
        seq = seq.replace('-','')
        print(strand)
        if strand == 1:
            aa_seq = Seq.translate(seq)
        if strand == -1:
            aa_seq = Seq.Seq(seq).reverse_complement().translate()
        else:    
            print('asdf')
        amino_acid_sequences[sample] = aa_seq

    return nucleotide_sequences,amino_acid_sequences

def get_genes_of_interest_for_all_patients(min_cov,min_base_freq):
    # builds meta data for multiple patients of itnerest

    genes_of_interest_dic = {}
    keys = patients
    value = None

    for patient in patients:
        pile_ups = load_pile_ups(patient)   
        coverage = get_coverage(pile_ups)  
        base_freq = get_base_frequency(pile_ups,coverage)   
        base_index = get_most_abundant_base(pile_ups)   
        diff_pos = get_different_positions(base_index,pile_ups) 
        filtred_pos = filter_diff_pos(pile_ups,coverage,base_freq,diff_pos,min_cov,min_base_freq)
        genes = get_genes(patient,filtred_pos)
        genes_of_interest_dic[patient] = genes
        

    return genes_of_interest_dic
    

def build_meta_single_patient(patient,min_cov,min_base_freq):
    # builds meta data for just one patient

    pile_ups = load_pile_ups(patient)   
    coverage = get_coverage(pile_ups)   
    base_freq = get_base_frequency(pile_ups,coverage)   
    base_index = get_most_abundant_base(pile_ups)   
    diff_pos = get_different_positions(base_index,pile_ups) 
    filtred_pos = filter_diff_pos(pile_ups,coverage,base_freq,diff_pos,min_cov,min_base_freq)
     
    return pile_ups, coverage, base_freq, base_index, diff_pos, filtred_pos
        
def store_positios_per_gene(genes_of_interest_dic):
    all_genes = []
    for patient in patients:
        for gene in genes_of_interest_dic[patient]:
            if gene[2] not in all_genes:
                all_genes.append(gene[2])

    pos_genes = {}            
    keys = all_genes
    values = None
    for key in keys:
        pos_genes[key] = values

    for patient in patients:
        pile_ups, coverage, base_freq, base_index, diff_pos, filtred_pos = build_meta_single_patient(patient,30,0.8)
        for gene in genes_of_interest_dic[patient]:
            nucleotide_sequences,amino_acid_sequences = translate_seq_from_pile_ups(patient,base_index,filtred_pos,gene)
            pos = a = [(pos,a,d) for pos,(a,d) in enumerate(zip(amino_acid_sequences['Sample0'], amino_acid_sequences['Sample1'])) if a!=d] #or a!=g or d!=g]
            pos_genes[gene[2]]=pos
    return pos_genes
                    
def rate_of_ESBL_genes(patient):
    pile_ups, coverage, base_freq, base_index, diff_pos, filtred_pos = build_meta_single_patient(patient,30,0.8)
    esbl_genes = ['ccdA_2','ccdB_2','finO','klcA_5','macB_3','resA','spo0C','traD','traI','xerD_6','ylpA']
    gb_file = base_dir+patient+'/Ref_genome/Annotations/prokka.gbk'
    for contig_counter,gb_record in enumerate(SeqIO.parse(open(gb_file,'r'),'genbank')):
        for esbl_gene in esbl_genes:
            for feature in gb_record.features:
                try:
                    if feature.qualifiers['gene'][0]==esbl_gene:
                        print('Ref_genome: Contig: ',contig_counter,'ESBL_gene: ',esbl_gene)
                except KeyError:
                    pass

    #check samples if they have the ESBL gene
    samples = sorted(os.listdir(base_dir+patient))[1:]
    for sample in samples:
        gb_file = base_dir+patient+'/'+sample+'/Spades/Annotations/prokka.gbk'
        print(gb_file)
        for contig_counter,gb_record in enumerate(SeqIO.parse(open(gb_file,'r'),'genbank')):
            for esbl_gene in esbl_genes:
                for feature in gb_record.features:
                    try:
                        if feature.qualifiers['gene'][0]==esbl_gene:
                            print('Contig: ',contig_counter,'Sample: ',sample,'ESBL_gene: ',esbl_gene)
                    except KeyError:
                        pass
    return pile_ups, coverage, base_freq, base_index, diff_pos, filtred_pos



    

    