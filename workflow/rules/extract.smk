# rule minimap2_index:
#     input:
#         target="target/{input1}.fasta"
#     output:
#         "{input1}.mmi"
#     log:
#         "logs/minimap2_index/{input1}.log"
#     params:
#         extra=""  # optional additional args
#     threads: 3
#     wrapper:
#         "0.66.0/bio/minimap2/index"
#
# rule minimap2:
#     input:
#         target="target/{input1}.mmi", # can be either genome index or genome fasta
#         query=["query/reads1.fasta", "query/reads2.fasta"]
#     output:
#         "aligned/{input1}_aln.paf"
#     log:
#         "logs/minimap2/{input1}.log"
#     params:
#         extra="-x map-pb"  # optional
#     threads: 3
#     wrapper:
#         "0.66.0/bio/minimap2/aligner"

rule edit_distance:  # calculates edit distances with minimap2
    input:
        "results/trimmed/{sample}.fastq.gz"
    output:
        "results/test/{sample}.test.txt"
    params:
        ""
    log:
        "logs/{sample}.test.log"
    conda:
        "../envs/mappy.yaml"
    script:
        "../scripts/edit_distance.py"

