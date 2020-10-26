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
        graph_xml="results/noderad/graph/{sample}.xml.gz",
        graph_figure="results/noderad/graph/{sample}.pdf"
    params:
        threshold=12, # threshold for maximum distance value for building the graph, default: 23
        mut_total=config["mutationrates"]["total"],
        mut_subst=config["mutationrates"]["substitution"],
        mut_ins=config["mutationrates"]["insertion"],
        mut_del=config["mutationrates"]["deletion"],
    log:
        "logs/noderad/graph/{sample}-graph.log"
    conda:
        "../envs/noderad_graph.yaml"
    script:
        "../scripts/noderad_graph.py"

# extracts the connected components and solves the ilp to determine the optimal representatives
rule noderad_representatives:
    input:
         "results/noderad/graph/{sample}.xml.gz"
    output:
        representatives="results/noderad/representatives/{sample}.txt",
        # optional output files:
        connected_components_xml="results/noderad/connected_components/{sample}.all_components.xml.gz",
        connected_components_figure="results/noderad/connected_components/{sample}.all_components.pdf",
        dir_subgraphs=directory("results/noderad/connected_components/subgraphes/{sample}")
    params:
        solver=config["ilp"]["solver"],
        mip=config["ilp"]["mip"],
        timelimit=config["ilp"]["timeLimit"],
        gaprel=config["ilp"]["gapRel"],
        gapabs=config["ilp"]["gapAbs"],
        fracgap=config["ilp"]["fracGap"],
        maxnodes=config["ilp"]["maxNodes"],
        maxmemory=config["ilp"]["maxMemory"]
    log:
        "logs/noderad/representatives/{sample}-representatives.log"
    threads:  # for the PuLP solver options the total number of threads must be divided among all samples
        int(config["ilp"]["threads"]/len(samples.index)) if config["ilp"]["threads"] else 1
    conda:
         "../envs/noderad_representatives.yaml"
    script:
        "../scripts/noderad_representatives.py"

