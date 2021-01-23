import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

sample = snakemake.wildcards.get('sample')

# input
vcf_in = pd.read_csv(snakemake.input.get("vcf_data"), sep="\t")[["#CHROM", "REF", "ALT"]]

individual = snakemake.params.get("individual").replace(" ", "_")

prefix = "{sample}|{indiv}|".format(sample=sample,indiv=individual)
fasta = ""

for i in range(len(vcf_in["#CHROM"])):
    name = vcf_in["#CHROM"][i]
    fasta += ">{prefix}{loc}-REF\n{ref}\n".format(prefix=prefix, loc=name, ref=vcf_in["REF"][i])
    variants = str(vcf_in["ALT"][i]).split(",")
    if not variants[0] == ".":
        idx = 1
        for var in variants:
            fasta += ">{prefix}{loc}-ALT_{idx}\n{alt}\n".format(prefix=prefix, loc=name, idx=idx, alt=var)
            idx += 1

# write fasta file
fasta_file = open(snakemake.output.get("vcf_fasta", ""), "w")
fasta_file.write(fasta)
