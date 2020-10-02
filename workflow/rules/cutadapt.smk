### data preprocessing ###

# demultiplexing and removing adapter sequences
rule cutadapt_se:
    input:
        config["resources"]["fastq-data"]
    output:
        fastq="results/trimmed/{sample}.fastq.gz",
        qc="results/trimmed/{sample}.se.qc.txt",
        log="logs/cutadapt/{sample}.se.log"#,
        # restfile="results/trimmed/restfiles/restfile-{sample}.fastq.gz"
    params:
        adapter=get_adapter,
        extra="--no-indels" # --action=lowercase -r {output.restfile}
    log:
        "logs/cutadapt/{sample}.se.log"
    conda:
        "../envs/temp_cutadapt_se.yaml"
    shell:
        # ToDo: move optional parameters into config.yaml
        "cutadapt {params.adapter} {params.extra} -o results/trimmed/{{name}}.fastq.gz {input} > {output.qc} 2> {log}"
    #TODO: possibly replace with wrapper after adaptation
    # wrapper:
    #     "xxxx/bio/cutadapt/se"
