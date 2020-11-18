rule minimap2_index:
    input:
        target="results/trimmed/{sample}.fastq.gz"
    output:
        "results/minimap2/{sample}.fastq.mmi"
    log:
        "logs/minimap2/{sample}.log"
    threads: 3
    wrapper:
        "0.66.0/bio/minimap2/index"

# determine edit distances and best alignment with minimap2
rule minimap2:
    input:
        target="results/minimap2/{sample}.fastq.mmi",
        query=["results/trimmed/{sample}.fastq.gz", "results/trimmed/{sample}.fastq.gz"]
    output:
        "results/minimap2/aligned/{sample}_aln.sam"
    log:
        "logs/minimap2/{sample}.log"
    params:
        extra="-a -A 10 --eqx --cs=long"
    threads: 3
    wrapper:
        "0.66.0/bio/minimap2/aligner"

# creates graph for nodeRAD calculations
rule noderad_graph:
    input:
        sam="results/minimap2/aligned/{sample}_aln.sam",
        fastq="results/trimmed/{sample}.fastq.gz"
    output:
        graph_xml="results/noderad/1_graph/{sample}.xml.gz",
        # optional output files:
        graph_figure="results/noderad/1_graph/{sample}.pdf"
    params:
        threshold_max_edit_distance=config["params"]["threshold_max_edit_distance"],
        mut_total=config["genome-properties"]["mutationrates"]["total"],
        mut_subst=config["genome-properties"]["mutationrates"]["substitution"],
        mut_ins=config["genome-properties"]["mutationrates"]["insertion"],
        mut_del=config["genome-properties"]["mutationrates"]["deletion"]
    log:
        "logs/noderad/1_graph/{sample}-graph.log"
    conda:
        "../envs/noderad_graph.yaml"
    script:
        "../scripts/noderad_graph.py"

# extracts the connected components and solves the ilp to determine the optimal representatives
rule noderad_allel_fractions:
    input:
         "results/noderad/1_graph/{sample}.xml.gz"
    output:
        connected_components_xml="results/noderad/2_allel_fractions/connected_components/{sample}.all_components.xml.gz",
        connected_components_figure="results/noderad/2_allel_fractions/connected_components/{sample}.all_components.pdf",
        dir_subgraphs=directory("results/noderad/2_allel_fractions/connected_components/subgraphes"+"/{sample}")
    params:
        # optional params
        ploidy=config["genome-properties"]["ploidy"],
        treshold_seq_noise=config["genome-properties"]["treshold-seq-noise"]
    log:
        "logs/noderad/2_allel_fractions/{sample}-representatives.log"
    conda:
         "../envs/noderad_allel_fractions.yaml"
    script:
        "../scripts/noderad_allel_fractions.py"


