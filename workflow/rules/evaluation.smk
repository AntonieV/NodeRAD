import os
from pathlib import Path

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
        "v1.25.0/bio/samtools/faidx"

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
        dbname= lambda wc, input: "{}/{}".format(os.path.dirname(input.fasta_sim[0]), Path(os.path.basename(input.fasta_sim[0])).stem)
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
