### quality control ###

rule fastqc:
    input:
        "results/trimmed/{sample}.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.66.0/bio/fastqc"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "results/qc/multiqc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.66.0/bio/multiqc"
