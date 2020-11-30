from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

####### DATA PREPROCESSING #######

####### load config and sample sheets #######

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype = str).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

report: "../report/workflow.rst"

####### wildcard constraints #######

wildcard_constraints:
    sample = "|".join(samples.index)

####### data preprocessing #######

def get_adapter(wildcards):
    return " -g {}={} ".format(wildcards.sample, samples.loc[wildcards.sample, "barcode_1"])

####### results #######
def get_multiqc_input(wildcards):
    multiqc_input = []
    for sample in samples.index:
        multiqc_input.extend (
            expand (
                [
                    "logs/cutadapt/{sample}.se.log",
                    "results/trimmed/{sample}.se.qc.txt"
                ],
                sample = sample
            )
        )

        multiqc_input.extend (
            expand (
                [
                    "results/qc/fastqc/{sample}_fastqc.zip",
                    "results/qc/fastqc/{sample}.html"
                ],
                sample = sample
            )
        )
    return multiqc_input

def all_input(wildcards):
    wanted_input = []

    # multiQC-report
    wanted_input.extend(["results/qc/multiqc/multiqc.html"])

    # trimming reads
    for sample in samples.index:
        # requires output: vcf file
        wanted_input.extend(expand(["results/noderad/4_vcf/{sample}.vcf"], sample=sample))

        # optional figures of graph-intermediate-steps:
        if config["graph-intermediate-steps"]["edit-distance-graph"]["activate-xml"]:
            wanted_input.extend(expand(["results/noderad/1_graph/{sample}.xml.gz"], sample=sample))

        if config["graph-intermediate-steps"]["edit-distance-graph"]["activate-figure"]:
            wanted_input.extend(expand(["results/noderad/1_graph/{sample}.pdf"], sample=sample))

        if config["graph-intermediate-steps"]["connected-components"]["graph"]["activate-xml"]:
            wanted_input.extend (
                expand (
                    [
                        "results/noderad/2_connected_components/{sample}.all_components.xml.gz"
                    ],
                    sample=sample
                )
            )

        if config["graph-intermediate-steps"]["connected-components"]["graph"]["activate-figure"]:
            wanted_input.extend (
                expand (
                    [
                        "results/noderad/2_connected_components/{sample}.all_components.pdf"
                    ],
                    sample=sample
                )
            )

        if config["graph-intermediate-steps"]["connected-components"]["subgraphs"]["activate"]:
            wanted_input.extend (
                directory(
                    expand (
                        [
                            "results/noderad/2_connected_components/subgraphes/{sample}"
                        ],
                        sample=sample
                    )
                )
            )

    return wanted_input
