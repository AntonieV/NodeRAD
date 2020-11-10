import sys
from graph_tool.all import *
import noderad_graph_operations
from scipy.stats import poisson
import numpy
import textwrap
import pysam


sys.stderr = open(snakemake.log[0], "w")


### input
graph = load_graph(snakemake.input.get("graph"))
alignment = snakemake.input.get("alignment")

### required output
repres_fasta = snakemake.output.get("repres_fasta", "")
clustered_reads = snakemake.output.get("clustered_reads", "")

### optional output
opt_repres_xml = snakemake.output.get("opt_repres_xml", "")
opt_repres_figure = snakemake.output.get("opt_repres_figure", "")
dir_clusters = snakemake.output.get("dir_clusters", "")

### set graph properties
for (name, prop) in noderad_graph_operations.get_properties("optimized-graph"):
    if name.startswith("v_"):
        vars()[name] = graph.vertex_properties[str(prop)]
    if name.startswith("e_"):
        vars()[name] = graph.edge_properties[str(prop)]

### set new graph properties
for (name, prop, prop_type) in noderad_graph_operations.get_new_properties("optimized-graph"):
    if name.startswith("v_"):
        vars()[name] = graph.new_vertex_property(prop_type)
        graph.vertex_properties[prop] = vars()[name]
    if name.startswith("e_"):
        vars()[name] = graph.new_edge_property(prop_type)
        graph.edge_properties[prop] = vars()[name]

n_edges = graph.num_edges()

### get connected components of the graph
all_components = noderad_graph_operations.get_components(graph)

representatives = []
for (comp, comp_nr) in all_components:
    nodes = list([comp.vertex_index[node] for node in comp.vertices()])
    max_prob = 0
    opt_repres = []
    for n in range(1, len(nodes) + 1):
        reprs = []
        probab_poisson_vals = []
        for node in comp.vertices():
            repr_n_vals = comp.vertex_properties['total-number-of-representatives'][node]
            if n in repr_n_vals:
                # poisson.pmf(k) = exp(-mu) * mu**k / k!
                # mu = d = sum-in-edges?, k = len(in_neighbors[node]? oder cov als gegeben?
                reprs.append(comp.vertex_index[node])
                predictand = comp.vertex_properties['sum-in-edges'][node]
                poisson_expected = poisson(predictand)
                d = len(comp.get_in_neighbors(comp.vertex_index[node]))
                probab_poisson_vals.append(poisson_expected.pmf(d))

        # product of Poisson expectation of sum of the incoming edges as preceding value
        if len(probab_poisson_vals) > 0:
            prob = numpy.prod(probab_poisson_vals)
            if prob > max_prob:
                max_prob = prob
                opt_repres = reprs

    # finds component vertex in graph and adds the graph vertex object to representatives list
    for r in opt_repres:
        node = find_vertex(graph, v_id, comp.vertex_properties["id"][r])[0]
        representatives.append(node)

### write fasta file and get ref_names and ref_length for header in sam file
for rep in representatives:
    # limits the sequence to 80 characters per line
    seq_partition = textwrap.fill(graph.vertex_properties["sequence"][rep], 80)
    with open(repres_fasta, 'a+') as f:
            print(">{}\n{}".format(graph.vertex_properties["name"][rep], seq_partition), file=f)

### get alignments of representatives from graph and write sam file, set properties for optimal solution
header_align = pysam.AlignmentFile(alignment, "rb")
with pysam.AlignmentFile(clustered_reads, "w", template=header_align) as f:
    reads = []
    for rep in representatives:
        v_sol_repr_neigh[rep] = 2
        for (query, ref) in graph.get_in_edges(rep):
            ed = graph.edge(query, ref)
            if v_sol_repr_neigh[query] == 0:
                v_sol_repr_neigh[query] = 1
            e_in_edges[ed] = True
            read = pysam.AlignedSegment()
            read.query_name = graph.edge_properties["sam-format-qname"][ed]
            read.query_sequence = graph.vertex_properties["sequence"][query]
            read.flag = graph.edge_properties["sam-format-flag"][ed]
            read.reference_id = int(graph.edge_properties["sam-format-rname"][ed])
            read.reference_start = graph.edge_properties["sam-format-pos"][ed]
            read.mapping_quality = graph.edge_properties["sam-format-mapq"][ed]
            read.cigar = list(eval(graph.edge_properties["sam-format-cigar"][ed]))
            read.next_reference_id = int(graph.edge_properties["sam-format-rnext"][ed])
            read.next_reference_start = graph.edge_properties["sam-format-pnext"][ed]
            read.template_length = graph.edge_properties["sam-format-tlen"][ed]
            qual = graph.vertex_properties["quality-q-vals"][query]
            if qual:
                read.query_qualities = qual
            read.tags = list(eval(graph.edge_properties["sam-format-all-tags"][ed]))
            reads.append(read)
    for read in reads:
        f.write(read)

### remove all edges that do not belong to the solution for graph figure
for node in graph.vertices():
    if graph.vertex_properties["is-representative-or-neighbor"][node] == 0:
        graph.clear_vertex(node)

sys.stderr.write("graph construction summary after removal of the edges that were not part of "
                 "optiomal solution:\n nodes:\t{}\n edges:\t{}\n "
                 "There were {} edges removed.\n "
                 "As optimal solution {} clusters were detected.".format(graph.num_vertices(),
                                                                                   graph.num_edges(),
                                                                                   n_edges - graph.num_edges(),
                                                                                   len(representatives)))

### save graph with optimal solution properties
if opt_repres_xml:
    graph.save(opt_repres_xml)

### draw graph
if opt_repres_figure:
    pos = sfdp_layout(graph)
    if opt_repres_figure:
        graph_draw(graph, vertex_fill_color=graph.vertex_properties['is-representative-or-neighbor'],
                   edge_color=graph.edge_properties['incoming-edges-to-representative'],
                   pos=pos, vertex_size=1, output=opt_repres_figure)

### draws each cluster separately
if dir_clusters:
    noderad_graph_operations.get_components(graph, wildcards=snakemake.wildcards.get('sample'),
                                            sub_dir=dir_clusters, v_color='is-representative-or-neighbor',
                                            e_color='incoming-edges-to-representative')

