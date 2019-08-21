# Patient isolates analysis

Sequencing data from the morbiodstat experiments 01 and 02 were analyzed with the scripts located in this folder. Sample 0 contains the sequencing data of the strains used for the experiments. 
Illumina sequencing data from morbidostat samples was first trimmed with submit_trim.sh and then analyzed for contaminations with bacillus_reads.py. Then the consensus sequence was hybrid assembled for the strains used for the morbidostat experiments with submit_unicycler.sh. Illumina data from the samples were mapped to the corresponding consensus sequences with mapping.sh and processed to pile ups with pileup.py. SNPs in morbidostat samples were identified with analysis_modular.py and annotated using the genbank files produced with prokka.sh.
