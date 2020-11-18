from graph_tool.all import *
import sys
import os
from itertools import combinations_with_replacement
from collections import Counter


def get_properties(graphtype):
    if graphtype == "distance-graph":
        return [(name, prop) for (name, prop, prop_type) in get_new_properties("init-graph")]


def get_new_properties(graphtype):
    if graphtype == "init-graph":
        return [("v_id", "id", "string"), ("v_name", "name", "string"),
                ("v_seq", "sequence", "string"),
                ("v_qual", "quality", "vector<float>"),  # p error of quality phred scores
                ("v_q_qual", "quality-q-vals", "vector<int>"),  # q values of quality phred scores
                ("e_dist", "distance", "int"),
                ("e_cs", "cs-tag", "string"),  # cigar string from cs tag, short or long option can be selected in minimap2 rule
                ("e_lh", "likelihood", "float")]

    if graphtype == "distance-graph":
        ### new properties for ILP calculation
        # v_concom: index number of the connected component to which the node belongs to
        # v_frac_lh: likelihood of allele fractions given one read
        # v_sum_in: sum of likelihood of all incoming edges of the node
        return [("v_concom", "component-label", "int32_t"),
                ("v_frac_lh", "read-fraction-likelihood", "float"),
                ("v_sum_in", "sum-in-edges", "double")]

def get_components(graph_obj, msg="", wildcards="", sub_dir=None, graph_xml=None, graph_figure=None, v_color=None, e_color=None):
    v_concom = graph_obj.vertex_properties["component-label"]
    all_components, hist = label_components(graph_obj, vprop=v_concom, directed=False)
    if msg:
        sys.stderr.write("\n\n{}:\n".format(msg))
        sys.stderr.write("Number of connected components: {}\n\n".format(max(all_components.a)))
        sys.stderr.write("Histogram of connected components:\n" + str(hist) + "\n\n")
        sys.stderr.write("connected components with more than one vertex are:\n\n")
    connected_components = []
    # default colors
    v_color_prop = "red"
    e_color_prop = "blue"

    for comp_nr in range(max(all_components)):
        # a is a shortcut for the get_array() method as an attribute,
        # please see https://graph-tool.skewed.de/static/doc/graph_tool.html#graph_tool.PropertyMap.get_array
        conn_component = GraphView(graph_obj, vfilt=all_components.a == comp_nr)
        conn_component = Graph(conn_component, prune=True)
        if hist[comp_nr] > 1:
            connected_components.append(tuple([conn_component, comp_nr]))
            if msg:
                sys.stderr.write("\t" + str(conn_component) + "\n")

            ### optional output files for each subgraph:
            if sub_dir:
                if not os.path.exists(sub_dir):
                    os.makedirs(sub_dir)
                conn_component.save(
                    "{}/{}-{}.subgraph.xml.gz".format(sub_dir, wildcards, str(comp_nr)))

                # set colors for nodes and edges and draw each component
                if v_color:
                    v_color_prop = conn_component.vertex_properties[v_color]
                if e_color:
                    e_color_prop = conn_component.edge_properties[e_color]

                pos = sfdp_layout(conn_component)
                graph_draw(conn_component, vertex_fill_color=v_color_prop, edge_color=e_color_prop, pos=pos,
                           vertex_size=9, output="{}/{}-{}.subgraph.pdf".format(sub_dir, wildcards, str(comp_nr)))

    ### optional output: all connected components of the graph:
    # set colors for nodes and edges and draw each component
    if v_color:
        v_color_prop = graph_obj.vertex_properties[v_color]
    if e_color:
        e_color_prop = graph_obj.edge_properties[e_color]
    if graph_xml:
        graph_obj.save(graph_xml)
    if graph_figure:
        pos = sfdp_layout(graph_obj)
        graph_draw(graph_obj, vertex_fill_color=v_color_prop,
                   edge_color=e_color_prop,
                   pos=pos, vertex_size=1, output=graph_figure)
    return connected_components


def get_candidate_alleles(graph, reads, noise):
    read_support = Counter()
    for read in reads:
        read_support[graph.vertex_properties['sequence'][read]] += 1
    return sorted(seq for seq, count in read_support.items() if count >= noise)

### returns a map with lists of all possible combinations of fractions from a given number of alleles n and
# the ploidy given in config file
def get_candidate_vafs(n, ploidy):
    n_alleles = n + n % ploidy  # each locus is of same ploidy
    get_combination_vafs = lambda comb: [sum(x == i for x in comb) / n_alleles for i in range(n)]
    return map(get_combination_vafs, combinations_with_replacement(range(n), n_alleles))


