import sys
import os
from graph_tool.all import *
from Bio import SeqIO
import likelihood_operations


def set_properties():
    return [("v_id", "id", "string"), ("v_name", "name", "string"),
            ("v_seq", "sequence", "string"),
            ("v_qual", "quality", "vector<float>"),  # p error of quality phred scores
            ("v_q_qual", "quality-q-vals", "vector<int>"),  # q values of quality phred scores
            ("v_concom", "component-label", "int32_t"),  # index number of the connected component to which the node belongs to
            ("e_dist", "distance", "int"),
            ("e_cs", "cs-tag", "string"),  # cigar string from cs tag, short or long option can be selected in minimap2 rule
            ("e_cigar_tup", "cigar-tuples", "string"),  # cigar tuples https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
            ("e_lh", "likelihood", "float"),
            ("v_concom", "component-label", "int32_t")]  # index number of the connected component to which the node belongs to


# add vertices
def set_nodes(_reads, graph, format, sample):
    for record in SeqIO.parse(_reads, format):
        node = graph.add_vertex()
        graph.vp["id"][node] = "{idx}_{sample}".format(idx=node, sample=sample)
        graph.vp["name"][node] = (record.id).split(' ', 1)[0]
        graph.vp["sequence"][node] = record.seq
        graph.vp["quality"][node] = [10 ** (-1 * i / 10) for i in record.letter_annotations["phred_quality"]]
        graph.vp["quality-q-vals"][node] = record.letter_annotations["phred_quality"]
    return graph


def set_edges(graph, read, threshold, stats, mut_rates):
    if read.has_tag("NM") and not read.query_name == read.reference_name:
        nm = 0  # edit distance between query and reference sequence
        cig = ""  # elements of cigar string
        for (tag, tag_val) in read.get_tags():
            if tag == "NM":
                nm = tag_val
            if tag == "cs":
                cig = tag_val
        if nm <= threshold:
            if stats["max_NM"] < nm:
                stats["max_NM"] = nm
            query_node = find_vertex(graph, graph.vp["name"], read.query_name)[0]
            ref_node = find_vertex(graph, graph.vp["name"], read.reference_name)[0]
            if not graph.vertex_index[ref_node] in graph.get_out_neighbors(query_node):  # ignores duplicates
                # add edge from alignment of sam-file
                edge = graph.add_edge(query_node, ref_node)
                graph.ep["distance"][edge] = nm
                graph.ep["cs-tag"][edge] = cig
                graph.ep["cigar-tuples"][edge] = read.cigar
                likelihood_calculation = likelihood_operations.get_alignment_likelihood(read.cigartuples,
                                                               graph.vertex_properties["quality"][query_node],
                                                               mut_rates, stats)
                graph.ep["likelihood"][edge] = likelihood_calculation[0]
                stats = likelihood_calculation[1]
    return graph, stats


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

            # optional output files for each subgraph:
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

    # optional output: all connected components of the graph:
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




