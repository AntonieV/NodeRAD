### quality control ###

rule fastqc:
    input:
        "results/trimmed/{sample}.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip"
    log:
        "logs/qc/fastqc/{sample}_fastqc.log"
    wrapper:
        "v1.25.0/bio/fastqc"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "results/qc/multiqc/multiqc.html"
    log:
        "logs/qc/multiqc.log"
    wrapper:
        "v1.25.0/bio/multiqc"
