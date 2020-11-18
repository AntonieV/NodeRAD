import sys
from graph_tool.all import *
import noderad_graph_operations


sys.stderr = open(snakemake.log[0], "w")


### input
graph = load_graph(snakemake.input[0])

### params
ploidy = snakemake.params.get("ploidy", "")
noise = snakemake.params.get("treshold_seq_noise", "")

if not noise:
    noise = 1

### optional output
connected_components_xml = snakemake.output.get("connected_components_xml", "")
connected_components_figure = snakemake.output.get("connected_components_figure", "")
dir_subgraphs = snakemake.output.get("dir_subgraphs", "")


### set graph properties
for (name, prop) in noderad_graph_operations.get_properties("distance-graph"):
    if name.startswith("v_"):
        vars()[name] = graph.vertex_properties[str(prop)]
    if name.startswith("e_"):
        vars()[name] = graph.edge_properties[str(prop)]

### set new graph properties
for (name, prop, prop_type) in noderad_graph_operations.get_new_properties("distance-graph"):
    if name.startswith("v_"):
        vars()[name] = graph.new_vertex_property(prop_type)
        graph.vertex_properties[prop] = vars()[name]
    if name.startswith("e_"):
        vars()[name] = graph.new_edge_property(prop_type)
        graph.edge_properties[prop] = vars()[name]


### pre-calculation of the sum of likelihoods of incoming edges for each node (to improve the runtime of the ILP)
for node in graph.vertices():
    likelihood_sum = 0
    for neighbor in graph.get_in_neighbors(graph.vertex_index[node]):
        likelihood = graph.edge_properties['likelihood'][graph.edge(neighbor, node)]
        likelihood_sum += likelihood
    v_sum_in[node] = likelihood_sum

### decoupling nodes whose sequences occur less frequently in the graph than the treshold_seq_noise
sys.stderr.write("Initial graph summary for sample {}:\n "
                 "nodes:\t{}\n edges:\t{}\n\n".format(snakemake.wildcards.get('sample'), graph.num_vertices(), graph.num_edges()))
filtered_seqs = noderad_graph_operations.get_candidate_alleles(graph, graph.vertices(), noise)
for node in graph.vertices():
    if not graph.vertex_properties['sequence'][node] in filtered_seqs:
        graph.clear_vertex(node)
sys.stderr.write("Graph summary for sample {} after filtering sequences that occur less than {} times :\n "
                 "nodes:\t{}\n edges:\t{}\n\n".format(snakemake.wildcards.get('sample'), noise, graph.num_vertices(), graph.num_edges()))


### extract connected components
message = "CONNECTED COMPONENTS based on the graph construction from the edit distances (minimap2)"
connected_components = noderad_graph_operations.get_components(graph, message, snakemake.wildcards.get('sample'),
                                                                   dir_subgraphs, connected_components_xml,
                                                                   connected_components_figure, v_color="component-label",
                                                                   e_color="likelihood")

### step 2
log_text = ""
for (comp, comp_nr) in connected_components:
    alleles = noderad_graph_operations.get_candidate_alleles(comp, comp.vertices(), noise)
    n = len(alleles)
    alleles_likelihood = dict(zip(alleles, [0] * n))
    fractions = list(noderad_graph_operations.get_candidate_vafs(n, ploidy))
    fraction_factors = set(sum(fractions, []))
    nodes = list([comp.vertex_index[node] for node in comp.vertices()])

    for node in nodes:
        fraction_likelihood_read = 0
        for factor in fraction_factors:
            fraction_likelihood_read += factor * comp.vertex_properties['sum-in-edges'][node]
        v_frac_lh[node] = fraction_likelihood_read
        alleles_likelihood[comp.vertex_properties['sequence'][node]] += fraction_likelihood_read

    sys.stderr.write("\nLikelihood of alleles in connected component {} in sample {}:\n\n".format(comp_nr, snakemake.wildcards.get('sample')))

    for allel in alleles_likelihood:
        sys.stderr.write("   likelihood = {} of sequence {}\n".format(alleles_likelihood[allel], allel))

    max = 0
    frac = []
    for fraction in fractions:
        fraction_likelihood = 1
        for i in range(n):
            if fraction[i] > 0:
                fraction_likelihood *= fraction[i] * alleles_likelihood[alleles[i]]
        if fraction_likelihood > max:
            max = fraction_likelihood
            frac = fraction
    log_text += "\n\nMaximal likelihood of allele fractions given all reads in connected component {} in sample {}:\n" \
                "   fraction likelihood = {}\n" \
                "   fraction: {}\n".format(comp_nr, snakemake.wildcards.get('sample'), max, frac)

sys.stderr.write(log_text)
