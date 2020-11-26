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

# required output
vcf = open(snakemake.output.get("vcf", ""), "w")
vcf_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(sample)
vcf.write(vcf_header)

# optional output
graph_xml = snakemake.output.get("graph_xml", "")
output_figure = snakemake.output.get("graph_figure", "")
connected_components_xml = snakemake.output.get("connected_components_xml", "")
connected_components_figure = snakemake.output.get("connected_components_figure", "")
dir_subgraphs = snakemake.output.get("components_subgraphs", "")
dir_alleles_subgraphs = snakemake.output.get("alleles_subgraphs", "")
dir_alleles_spanning_trees = snakemake.output.get("alleles_spanning_trees")

# params for mutation rates, ploidy, threshold and noise
threshold = snakemake.params.get("threshold_max_edit_distance", "")
ploidy = snakemake.params.get("ploidy", "")

noise = snakemake.params.get("treshold_seq_noise", "")

# get reads format
format = ""  # format of reads
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
for (name, prop, prop_type) in graph_operations.set_properties("reads-alignment"):
    if name.startswith("g_"):
        vars()[name] = graph.new_graph_property(prop_type)
        graph.graph_properties[prop] = vars()[name]
    if name.startswith("v_"):
        vars()[name] = graph.new_vertex_property(prop_type)
        graph.vertex_properties[prop] = vars()[name]
    if name.startswith("e_"):
        vars()[name] = graph.new_edge_property(prop_type)
        graph.edge_properties[prop] = vars()[name]

# set mutation-rates and heterozygosity
graph.graph_properties["subst-mutation-rates"] = snakemake.params.get("mut_subst", "")
graph.graph_properties["ins-mutation-rates"] = snakemake.params.get("mut_ins", "")
graph.graph_properties["del-mutation-rates"] = snakemake.params.get("mut_del", "")

graph.graph_properties["subst-heterozygosity"] = snakemake.params.get("heterozyg_subst", "")
graph.graph_properties["ins-heterozygosity"] = snakemake.params.get("heterozyg_ins", "")
graph.graph_properties["del-heterozygosity"] = snakemake.params.get("heterozyg_del", "")

# create first empty node for graph
node = graph.add_vertex()
v_id[node] = "{idx}_{sample}".format(idx=0, sample=sample)
v_name[node] = ""
v_seq[node] = ""
v_qual[node] = ""

# add reads as vertices of the graph
if reads.endswith(".gz"):
    with gzip.open(reads, "rt") as _reads:
        graph = graph_operations.set_nodes(graph, _reads, format, sample)
else:
    with open(reads, "rU") as _reads:
        graph = graph_operations.set_nodes(graph, _reads, format, sample)

# add edges from all-vs-all alignment of reads (please see rule minimap2)
stats = {"max_NM": 0, "n_snp": 0, "n_ins": 0, "n_del": 0}
sam = pysam.AlignmentFile(sam, "rb")
for read in sam.fetch(until_eof=True):
    graph_build = graph_operations.set_edges(graph, read, threshold, stats)
    graph = graph_build[0]
    stats = graph_build[1]
sam.close()
graph.remove_vertex(0)

# statistics
max_NM = stats["max_NM"]  # highest edit distance
n_snp = stats["n_snp"]  # counter for total number of substitutions/snp's
n_ins = stats["n_ins"]  # counter for total number of insertions
n_del = stats["n_del"]  # counter for total number of deletions

# write log files
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

graph_operations.save_and_draw_graph(graph, xml_out=graph_xml,
                                     figure_out=output_figure,
                                     e_color="likelihood")

# in case these directories are expected as optional output, but some samples do not produce
# subgraphs, e.g. spanning trees. This way the empty directories for these samples preserved
# and are not removed by snakemake.
graph_operations.set_dirs([dir_subgraphs, dir_alleles_subgraphs, dir_alleles_spanning_trees])

# extract connected components
message = "CONNECTED COMPONENTS based on the graph construction from the edit distances (minimap2)"
connected_components = graph_operations.get_components(graph, message, snakemake.wildcards.get('sample'),
                                                       dir_subgraphs, connected_components_xml,
                                                       connected_components_figure, v_color="component-label",
                                                       e_color="likelihood")

loc_nr = 1
for (comp, comp_nr) in connected_components:
    # step 2: likelihood of allele fractions
    alleles = likelihood_operations.get_candidate_alleles(comp, comp.vertices(), noise)
    n = len(alleles)
    vafs_candidates = list(likelihood_operations.get_candidate_vafs(n, ploidy))
    nodes = list([comp.vertex_index[node] for node in comp.vertices()])

    # calculate the likelihood over ALL reads
    vafs_likelihoods = [likelihood_operations.calc_vafs_likelihood(comp, vafs, nodes, alleles) for vafs in
                        vafs_candidates]
    if not vafs_likelihoods:  # case empty list
        continue
    max_likelihood_idx = np.argmax(vafs_likelihoods)

    # obtain ML solution for step 2
    max_likelihood_vafs = vafs_candidates[max_likelihood_idx]

    # write to log file
    sys.stderr.write(
        "\n\nStats for component {} in sample {} with {} alleles and ploidy = {}:\n".format(comp_nr, sample,
                                                                                              n, ploidy))
    sys.stderr.write("\n\tMaximum vafs likelihood:\n")
    for vaf, allele in zip(max_likelihood_vafs, alleles):
        sys.stderr.write("\t\t{} for allele {}\n".format(vaf, allele))

    # step 3: likelihood of loci given alleles and fractions
    loci_candidates = list(likelihood_operations.get_candidate_loci(n, ploidy))
    allele_subgraph = likelihood_operations.get_allele_subgraph(comp, alleles, dir_alleles_subgraphs, comp_nr)
    spanning_tree = likelihood_operations.get_spanning_tree(allele_subgraph, dir_alleles_spanning_trees, comp_nr)

    loci_likelihoods = [
        likelihood_operations.calc_loci_likelihoods(comp, max_likelihood_vafs, alleles, loci, allele_subgraph,
                                                    spanning_tree)
        for loci in loci_candidates
    ]
    # write to log file
    sys.stderr.write("\n\tMaximum loci likelihood calculations:\n")
    max_likelihood_idx = np.argmax(loci_likelihoods)
    max_likelihood_loci = loci_candidates[max_likelihood_idx]
    sys.stderr.write("\t\tloci_likelihoods:\n\t\t\t{}\n\t\tmax_likelihood_idx:\n\t\t\t{}"
                     "\n\t\tmax_likelihood_loci:\n\t\t\t{}\n".format(loci_likelihoods, max_likelihood_idx,
                                                                    max_likelihood_loci))

    # step 4: results output to VCF file
    loci_alleles = likelihood_operations.get_sorted_loci_alleles(alleles, max_likelihood_loci)
    n_alleles = len(loci_alleles)
    genotype = likelihood_operations.get_genotype(n_alleles)
    alt = ', '.join(loci_alleles[1:]) if n_alleles > 1 else "."
    if n_alleles:
        vcf.write("{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tGT\t{gt}\n".format(chrom="LOC{}".format(loc_nr),
                                                                             pos="1", ref=loci_alleles[0],
                                                                             alt=alt, gt=genotype))
        loc_nr += 1
