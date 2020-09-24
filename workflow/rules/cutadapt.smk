### data preprocessing ###

# demultiplexing and removing adapter sequences
rule cutadapt_se:
    input:
        config["resources"]["fastq-data"]
    output:
        fastq="results/trimmed/{sample}.fastq.gz",
        qc="results/trimmed/{sample}.se.qc.txt",
        log="logs/cutadapt/{sample}.se.log",
        restfile="results/trimmed/restfiles/restfile-{sample}.fastq.gz"
    params:
        get_adapter
    log:
        "logs/cutadapt/{sample}.se.log"
    conda:
        "../envs/temp_cutadapt_se.yaml"
    shell:
        # ToDo: move optional parameters into config.yaml
        # --action=lowercase -r {output.restfile}
        "cutadapt {params} --no-indels -o results/trimmed/{{name}}.fastq.gz {input} > {output.qc} 2> {log}"
    #TODO: possibly replace with wrapper after adaptation
    # wrapper:
    #     "xxxx/bio/cutadapt/se"
