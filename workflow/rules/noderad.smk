rule minimap2_index:
    input:
        target="results/trimmed/{sample}.fastq.gz"
    output:
        "results/minimap2/{sample}_aln.mmi"
    log:
        "logs/minimap2/{sample}.log"
    threads: 3
    wrapper:
        "0.66.0/bio/minimap2/index"

# determine edit distances and best alignment with minimap2
rule minimap2:
    input:
        target="results/minimap2/{sample}_aln.mmi",
        query=["results/trimmed/{sample}.fastq.gz", "results/trimmed/{sample}.fastq.gz"]
    output:
        pipe("results/minimap2/{sample}_aln.sam")
    log:
        "logs/minimap2/{sample}.log"
    params:
        extra="-a -A 10 -N 100 --eqx --cs=long"
    threads: 3
    wrapper:
        "0.66.0/bio/minimap2/aligner"

# converts SAM files to BAM files
rule samtools_view:
    input:
        "results/minimap2/{sample}_aln.sam"
    output:
        "results/minimap2/{sample}_aln.bam"
    params:
        " -b "
    wrapper:
        "v0.69.0/bio/samtools/view"

# RADSeq analysis: calulates alleles and loci likelihoods and returns vcf file with with the most probable loci
rule noderad:
    input:
        bam="results/minimap2/{sample}_aln.bam",
        fastq="results/trimmed/{sample}.fastq.gz"
    output:
        vcf="results/noderad/4_vcf/{sample}.vcf",
        # optional output files:
        graph_xml="results/noderad/1_graph/{sample}.xml.gz",
        graph_figure="results/noderad/1_graph/{sample}.pdf",
        connected_components_xml="results/noderad/2_connected_components/{sample}.all_components.xml.gz",
        connected_components_figure="results/noderad/2_connected_components/{sample}.all_components.pdf",
        components_subgraphs=directory("results/noderad/2_connected_components/subgraphes/{sample}")
    params:
        threshold_max_edit_distance=config["params"]["threshold_max_edit_distance"],
        # a diploid chromosome set is determined for this prototype,
        # for future use it can be configured in config["genome-properties"]["ploidy"]
        ploidy=config["genome-properties"]["ploidy"],
        treshold_seq_noise_small=config["genome-properties"]["noise"]["treshold-seq-noise"]["small-clusters"],
        treshold_seq_noise_large=config["genome-properties"]["noise"]["treshold-seq-noise"]["large-clusters"],
        treshold_cluster_size=config["genome-properties"]["noise"]["treshold-cluster-size"],
        # sequencing errors
        err_ins=config["genome-properties"]["errors-per-base"]["insertion"],
        err_del=config["genome-properties"]["errors-per-base"]["deletion"],
        # heterozygosity
        heterozyg_subst=config["genome-properties"]["heterozygosity"]["substitution"],
        heterozyg_ins=config["genome-properties"]["heterozygosity"]["insertion"],
        heterozyg_del=config["genome-properties"]["heterozygosity"]["deletion"],

    log:
        "logs/noderad/{sample}.log"
    conda:
        "../envs/noderad.yaml"
    script:
        "../scripts/noderad_main.py"

# evaluation of results

rule simulated_data_to_fasta:
    input:
         sim_data_stats=config["eval-data"]["small-dataset"]
    output:
        fasta_sim="results/evaluation/sim_fasta/{sample}.sim.fasta"
    params:
        individual=lambda wc: get_individual(wc.sample)
    log:
        "logs/evaluation/{sample}.sim_to_fasta.log"
    conda:
        "../envs/sim_to_fasta.yaml"
    script:
        "../scripts/sim_to_fasta.py"

