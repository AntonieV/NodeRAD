import os
from pathlib import Path

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
        heterozyg_del=config["genome-properties"]["heterozygosity"]["deletion"]
    log:
        "logs/noderad/{sample}.log"
    conda:
        "../envs/noderad.yaml"
    script:
        "../scripts/noderad_main.py"

# evaluation of results
rule simulated_data_to_fasta:
    input:
        sim_data_stats=config["eval-data"]
    output:
        fasta_sim="results/evaluation/sim_fasta/{sample}/{sample}.sim.fasta"
    params:
        individual=lambda wc: get_individual(wc.sample)
    log:
        "logs/evaluation/parse/{sample}.sim_to_fasta.log"
    conda:
        "../envs/sim_to_fasta.yaml"
    script:
        "../scripts/sim_to_fasta.py"

rule vcf_to_fasta:
    input:
        vcf_data="results/noderad/4_vcf/{sample}.vcf"
    output:
        vcf_fasta="results/evaluation/vcf_fasta/{sample}/{sample}.vcf.fasta"
    params:
        individual=lambda wc: get_individual(wc.sample)
    log:
        "logs/evaluation/parse/{sample}.vcf_to_fasta.log"
    script:
        "../scripts/vcf_to_fasta.py"

rule samtools_index:
    input:
        expand("results/evaluation/{type}_fasta/{sample}/{sample}.{type}.fasta", sample=samples.index, type=["sim", "vcf"])
    output:
        "results/evaluation/{type}_fasta/{sample}/{sample}.{type}.fasta.fai"
    params:
        ""
    log:
        "logs/evaluation/index/{type}_fasta/{sample}.{type}.fasta"
    wrapper:
        "v0.69.0/bio/samtools/faidx"

rule blast_database:
    input:
        fasta_sim="results/evaluation/sim_fasta/{sample}/{sample}.sim.fasta",
        fasta_idx="results/evaluation/sim_fasta/{sample}/{sample}.sim.fasta.fai"
    output:
        multiext("results/evaluation/sim_fasta/{sample}/{sample}.sim.fasta.",
                 "ndb", "nos", "not", "ntf", "nto")
    params:
        ""
    log:
        "logs/evaluation/blast/{sample}.blast_db.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "makeblastdb -in {input.fasta_sim} -input_type fasta -dbtype nucl -parse_seqids"

rule blast:
    input:
        fasta_sim=multiext("results/evaluation/sim_fasta/{sample}/{sample}.sim.fasta.",
                 "ndb", "nos", "not", "ntf", "nto"),
        res="results/evaluation/vcf_fasta/{sample}/{sample}.vcf.fasta"
    output:
        blast=pipe("results/evaluation/blast/{sample}.blast.raw")
    params:
        percent_identity=80,
        dbname= lambda wc, input: os.path.dirname(input.fasta_sim[0])+"/"+Path(os.path.basename(input.fasta_sim[0])).stem
    log:
        "logs/evaluation/blast/{sample}.log"
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -query {input.res} -outfmt 6 -db {params.dbname} -out {output.blast}"

rule blast_header:
    input:
        "results/evaluation/blast/{sample}.blast.raw"
    output:
        blast="results/evaluation/blast/{sample}.blast.tsv"
    params:
        header="query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\n"
    log:
        "logs/evaluation/blast/{sample}.header.log"
    shell:
        "echo -ne \"{params.header}\" | cat - {input} > {output}"

rule plots_blast:
    input:
        "results/evaluation/blast/{sample}.blast.tsv"
    output:
        ident=report("results/evaluation/blast/plots/perc_ident/{sample}.plot_loci.pdf", caption="../report/plot_perc_ident.rst", category="Blast"),
        ident_hist=report("results/evaluation/blast/plots/hist_perc_ident/{sample}.plot_hist.pdf", caption="../report/plot_hist_perc_ident.rst", category="Blast"),
        bit_scores=report("results/evaluation/blast/plots/bitscores/{sample}.plot_bitscores.pdf", caption="../report/plot_bitscores.rst", category="Blast"),
        evalues=report("results/evaluation/blast/plots/evalues/{sample}.plot_evalues.pdf", caption="../report/plot_evalues.rst", category="Blast")
    params:
        ""
    log:
        "logs/evaluation/blast/{sample}.plot_ident.log"
    conda:
        "../envs/plots_blast.yaml"
    script:
        "../scripts/plots_blast.R"
