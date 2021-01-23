import sys
import yaml
from itertools import chain
import re

sys.stderr = open(snakemake.log[0], "w")

sample = snakemake.wildcards.get('sample')

# input
sim_in = snakemake.input.get("sim_data_stats")
individual_name=snakemake.params.get("individual")
individual = snakemake.params.get("individual").replace("_", " ")

fasta = ""

# read yaml file
yaml_in = open(sim_in, 'r')
yaml_in_load = list(yaml.load_all(yaml_in, Loader=yaml.FullLoader))
sim_data_file = yaml_in_load[1]

# parse fasta data from yaml
for key_loc in sim_data_file.keys():
    # loci sequences
    seq = sim_data_file[key_loc]["p5 seq"]

    # loci mutations of given individual
    indiv_data = sim_data_file[key_loc]["individuals"][individual]
    mutations = list(chain(*[indiv_data[item]["mutations"] for item in indiv_data if indiv_data[item]["mutations"]]))

    # modifying the sequence when mutations are present for given individual and locus
    if mutations:
        mut_list = [re.split("@", item)[1] for item in mutations]
        mut_tuples = [tuple(re.split(":", mut)) for mut in mut_list]
        for (pos, mut) in mut_tuples:
            idx = pos
            if '(' or ')' in pos:
                idx = re.split("\(", re.split("\)", pos)[0])[-1]
            idx = int(idx) - 1
            if ">" in mut:
                mut_subs = re.split(">", mut)
                if seq[idx] == mut_subs[0]:
                    seq = seq[:idx] + mut_subs[-1] + seq[idx + 1:]
            if "+" in mut:
                mut_subs = re.split('\+', mut)[-1]
                seq = seq[:idx] + mut_subs + seq[idx:]
            if "-" in mut:
                mut_subs = re.split('\-', mut)[-1]
                if all([seq[idx] == mut_subs[i] for i in range(len(mut_subs))]):
                    seq = seq[:idx] + seq[idx + len(mut_subs):]

    # parsing to fasta format
    fasta_id = ">{}|{}|simulated".format(individual_name, key_loc.replace(" ", "_"))
    fasta += "{}\n{}\n".format(fasta_id, seq)

# write fasta file
fasta_file = open(snakemake.output.get("fasta_sim", ""), "w")
fasta_file.write(fasta)
