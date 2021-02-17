# The main entry point of workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

# The data preprocessing is partially based on the following implementations:
#        * https://github.com/snakemake-workflows/dna-seq-varlociraptor
#        * https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth
#        * https://github.com/snakemake-workflows/chipseq

# The data preprocessing steps include (cutadapt.smk and qc.smk):
#        * simulation of data for testing was performed with ddRAGE,
#          for more information please see https://ddrage.readthedocs.io/en/latest/getting-started/
#        * raw read quality control with FastQC and MultiQC
#        * adapter trimming with cutadapt

# The data analysis steps with NodeRAD include (noderad.smk):
#        * determine edit distances with minimap2 (rules minimap2_index and minimap2)
#        * convert output to .bam format with samtools view (rule samtools view)
#        * construct graph from edit distances with graph-tool (rule noderad)
#        * get connected components of the graph with graph-tool (rule noderad)
#        * calculation of the allele fractions and their likelihood considering the sequencing error rate (rule noderad)
#        * calculation of the most probable loci assignment from the allele fraction with the highest likelihood considering
#          the probabilities of heterozygosity (rule noderad)
#        * write vcf file with the identified loci of each connected component (rule noderad)

# The data evaluation steps include (evaluation.smk):
#        * parse the information (.yaml) of simulated data set from ddRAGE to FASTA format (rule simulated_data_to_fasta)
#        * parse the vcf file from NodeRAD analysis to FASTA format (rule vcf_to_fasta)
#        * indexing FASTA file from ddRAGE with samtools index (rule samtools_index) and create a BLAST database (rule blast_database)
#        * BLAST analysis of the loci identified by NodeRAD versus the loci from the simulation with ddRAGE (rule blast)
#        * add header to results of BLAST analysis (rule blast_header) and create plots (rule plots_blast)


include: "workflow/rules/common.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/cutadapt.smk"
include: "workflow/rules/noderad.smk"
include: "workflow/rules/evaluation.smk"

rule all:
    input: all_input
