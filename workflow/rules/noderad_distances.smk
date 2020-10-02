rule minimap2_index:
    input:
        target="results/trimmed/{sample}.fastq.gz"
    output:
        "results/minimap2/{sample}.fastq.mmi"
    log:
        "logs/minimap2_index/{sample}.log"
    threads: 3
    wrapper:
        "0.66.0/bio/minimap2/index"

rule minimap2:
    input:
        target="results/minimap2/{sample}.fastq.mmi",
        query=["results/trimmed/{sample}.fastq.gz", "results/trimmed/{sample}.fastq.gz"]
    output:
        "results/minimap2/aligned/{sample}_aln.sam"
    log:
        "logs/minimap2/{sample}.log"
    params:
        extra="-aY --eqx -B 0 --cs=long"
    threads: 3
    wrapper:
        "0.66.0/bio/minimap2/aligner"

# creates graph for nodeRAD calculations
rule create_noderad_graph:
    input:
        sam="results/minimap2/aligned/{sample}_aln.sam",
        fastq="results/trimmed/{sample}.fastq.gz"
    output:
        graph_xml="results/noderad_graph/{sample}.xml.gz",
        graph_figure="results/noderad_graph/{sample}.pdf"
    params:
        threshold=12 # threshold for maximum distance value for building the graph, default: 23
    log:
        "logs/noderad_graph/{sample}.test.log"
    conda:
        "../envs/noderad_graph.yaml"
    script:
        "../scripts/noderad_graph.py"

