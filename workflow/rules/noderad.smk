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
         "results/noderad/1_graph/{sample}.xml.gz"
    output:
        representatives="results/noderad/3_representatives/{sample}.tsv",
        representatives_xml="results/noderad/3_representatives/{sample}.xml.gz",
        # optional output files:
        connected_components_xml="results/noderad/2_connected_components/{sample}.all_components.xml.gz",
        connected_components_figure="results/noderad/2_connected_components/{sample}.all_components.pdf",
        dir_subgraphs=directory("results/noderad/2_connected_components/subgraphes"+"/{sample}"),
        repesentatives_figure="results/noderad/3_representatives/{sample}.pdf"
    params:
        # optional params for PuLP solvers
        solver=config["ilp"]["solver"],
        mip=config["ilp"]["mip"],
        timelimit=config["ilp"]["timeLimit"],
        gaprel=config["ilp"]["gapRel"],
        gapabs=config["ilp"]["gapAbs"],
        maxnodes=config["ilp"]["maxNodes"],
        maxmemory=config["ilp"]["maxMemory"]
    log:
        "logs/noderad/representatives/{sample}-representatives.log"
    threads:
        # for PuLP solver params the total number of threads must be divided among all samples
        int(config["ilp"]["threads"]/len(samples.index)) if config["ilp"]["threads"] else 1
    conda:
         "../envs/noderad_representatives.yaml"
    script:
        "../scripts/noderad_representatives.py"

# determine the best solution by considering distribution and coverage
rule optimized_solution:
    input:
        graph="results/noderad/3_representatives/{sample}.xml.gz",
        alignment="results/minimap2/aligned/{sample}_aln.sam"
    output:
        repres_fasta="results/noderad/4_optimized_solution/{sample}.fasta",
        clustered_reads="results/noderad/4_optimized_solution/{sample}.sam",
        # optional output files
        opt_repres_xml="results/noderad/4_optimized_solution/{sample}.xml.gz",
        opt_repres_figure="results/noderad/4_optimized_solution/{sample}.pdf",
        dir_clusters=directory("results/noderad/4_optimized_solution/cluster"+"/{sample}")
    params:
        ""
    log:
        "logs/noderad/optimized_solution/{sample}-representatives.log"
    conda:
         "../envs/noderad_optimized_solution.yaml"
    script:
        "../scripts/noderad_optimized_solution.py"
