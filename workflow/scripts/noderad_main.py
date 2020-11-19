import sys
import gzip
import pysam
from graph_tool.all import *
import numpy as np
import graph_operations
import likelihood_operations



sys.stderr = open(snakemake.log[0], "w")

sample = snakemake.wildcards.get('sample')

# input
sam = snakemake.input.get("sam")
reads = snakemake.input.get("fastq")

# optional output
output_figure = snakemake.output.get("graph_figure", "")
connected_components_xml = snakemake.output.get("connected_components_xml", "")
connected_components_figure = snakemake.output.get("connected_components_figure", "")
dir_subgraphs = snakemake.output.get("dir_subgraphs", "")

# params for mutation rates, ploidy, threshold and noise
threshold = snakemake.params.get("threshold_max_edit_distance", "")
ploidy = snakemake.params.get("ploidy", "")
mut_rates = {"mut_total": snakemake.params.get("mut_total", ""), "mut_subst": snakemake.params.get("mut_subst", ""),
             "mut_ins": snakemake.params.get("mut_ins", ""), "mut_del": snakemake.params.get("mut_del", "")}
noise = snakemake.params.get("treshold_seq_noise", "")

# get reads format
format = ""    # format of reads
if reads.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
    format = "fastq"
if reads.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")):
    format = "fasta"

# default value of threshold for maximum distance value
if threshold == "":
    threshold = 23
else:
    threshold = int(threshold)

# default value of noise
if not noise:
    noise = 1

# init graph
graph = Graph(directed=True)

# set graph properties
for (name, prop, prop_type) in graph_operations.set_properties():
    if name.startswith("v_"):
        vars()[name] = graph.new_vertex_property(prop_type)
        graph.vertex_properties[prop] = vars()[name]
    if name.startswith("e_"):
        vars()[name] = graph.new_edge_property(prop_type)
        graph.edge_properties[prop] = vars()[name]

# create first empty node for graph
node = graph.add_vertex()
v_id[node] = "{idx}_{sample}".format(idx=0, sample=sample)
v_name[node] = ""
v_seq[node] = ""
v_qual[node] = ""

# add reads as vertices of the graph
if reads.endswith(".gz"):
    with gzip.open(reads, "rt") as _reads:
        graph = graph_operations.set_nodes(_reads, graph, format, sample)
else:
    with open(reads, "rU") as _reads:
        graph = graph_operations.set_nodes(_reads, graph, format, sample)

# add edges from all-vs-all alignment of reads (please see rule minimap2)
stats = {"max_NM": 0, "n_snp": 0, "n_ins": 0, "n_del": 0}
sam = pysam.AlignmentFile(sam, "rb")
for read in sam.fetch(until_eof=True):
    graph_build = graph_operations.set_edges(graph, read, threshold, stats, mut_rates)
    graph = graph_build[0]
    stats = graph_build[1]
sam.close()
graph.remove_vertex(0)

# statistics
max_NM = stats["max_NM"]  # highest edit distance
n_snp = stats["n_snp"]   # counter for total number of substitutions/snp's
n_ins = stats["n_ins"]   # counter for total number of insertions
n_del = stats["n_del"]   # counter for total number of deletions

sys.stderr.write(
    "graph construction summary for sample {}:"
    "\n nodes:\t{}\n edges:\t{}\n max distance:\t{}".format(sample,
                                                            graph.num_vertices(),
                                                            graph.num_edges(),
                                                            stats["max_NM"]))
sys.stderr.write(
    "\n\tmutations:\n\t\tsnp's and substitutions: {}"
    "\n\t\tinsertions: {}\n\t\tdeletions: {}\n\n".format(stats["n_snp"],
                                                         stats["n_ins"],
                                                         stats["n_del"]))

graph.save(snakemake.output.get("graph_xml"))
pos = sfdp_layout(graph)
if output_figure:
    graph_draw(graph, vertex_color=[1, 1, 1, 0],
               edge_color=graph.edge_properties["likelihood"],
               pos=pos, vertex_size=1, output=output_figure)


# extract connected components
message = "CONNECTED COMPONENTS based on the graph construction from the edit distances (minimap2)"
connected_components = graph_operations.get_components(graph, message, snakemake.wildcards.get('sample'),
                                                                   dir_subgraphs, connected_components_xml,
                                                                   connected_components_figure, v_color="component-label",
                                                                   e_color="likelihood")

# step 2
for (comp, comp_nr) in connected_components:
    alleles = likelihood_operations.get_candidate_alleles(comp, comp.vertices(), noise)
    n = len(alleles)
    vafs_candidates = list(likelihood_operations.get_candidate_vafs(n, ploidy))
    nodes = list([comp.vertex_index[node] for node in comp.vertices()])

    # calculate the likelihood over ALL reads
    vafs_likelihoods = [likelihood_operations.calc_vafs_likelihood(vafs, comp, nodes, alleles, mut_rates) for vafs in vafs_candidates]
    max_likelihood_idx = np.argmax(vafs_likelihoods)
    # obtain ML solution for step 2
    max_likelihood_vafs = vafs_candidates[max_likelihood_idx]
    print(max_likelihood_vafs)

    sys.stderr.write("\nvafs with maximum likelihood in sample {}, copmonent {} is:\n".format(sample, comp_nr))
    for vaf, allele in zip(max_likelihood_vafs, alleles):
        sys.stderr.write("    {} for allele {}\n".format(vaf, allele))


    # step 3
    loci_candidates = likelihood_operations.get_candidate_loci(n, ploidy)

    # loci_likelihoods = [
    #     likelihood_operations.calc_loci_likelihoods(max_likelihood_vafs, comp, alleles, loci) for loci in loci_candidates
    # ]
    # max_likelihood_idx = np.argmax(loci_likelihoods)
    # max_likelihood_loci = loci_candidates[max_likelihood_idx]

