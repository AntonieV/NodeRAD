import sys
import os
from graph_tool.all import *
from Bio import SeqIO
import likelihood_operations


def set_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def set_properties(graph_type):
    if graph_type == "reads-alignment":
        return [("g_total_mut", "total-mutation-rates", "float"),
                ("g_subst_mut", "subst-mutation-rates", "float"),
                ("g_ins_mut", "ins-mutation-rates", "float"),
                ("g_del_mut", "del-mutation-rates", "float"),
                ("g_total_heterozyg", "total-heterozygosity", "float"),
                ("g_subst_heterozyg", "subst-heterozygosity", "float"),
                ("g_ins_heterozyg", "ins-heterozygosity", "float"),
                ("g_del_heterozyg", "del-heterozygosity", "float"),
                ("v_id", "id", "string"), ("v_name", "name", "string"),
                ("v_seq", "sequence", "string"),
                ("v_qual", "quality", "vector<float>"),  # p error of quality phred scores
                ("v_q_qual", "quality-q-vals", "vector<int>"),  # q values of quality phred scores
                ("v_concom", "component-label", "int32_t"),  # index number of the connected component to which the node belongs to
                ("e_dist", "distance", "int"),
                ("e_cs", "cs-tag", "string"),  # cigar string from cs tag, short or long option can be selected in minimap2 rule
                ("e_cigar_tup", "cigar-tuples", "string"),  # cigar tuples https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
                ("e_lh", "likelihood", "float")]
    if graph_type == "alleles":
        return [("g_total_heterozyg", "total-heterozygosity", "float"),
                ("g_subst_heterozyg", "subst-heterozygosity", "float"),
                ("g_ins_heterozyg", "ins-heterozygosity", "float"),
                ("g_del_heterozyg", "del-heterozygosity", "float"),
                ("v_al_seq", "allele-sequence", "string"),
                ("e_al_dist", "allele-distance", "int"),
                ("e_al_cigar_tup", "allele-cigar-tuples", "string")]


# add vertices
def set_nodes(graph, _reads, format, sample):
    for record in SeqIO.parse(_reads, format):
        node = graph.add_vertex()
        graph.vp["id"][node] = "{idx}_{sample}".format(idx=node, sample=sample)
        graph.vp["name"][node] = (record.id).split(' ', 1)[0]
        graph.vp["sequence"][node] = record.seq
        graph.vp["quality"][node] = [10 ** (-1 * i / 10) for i in record.letter_annotations["phred_quality"]]
        graph.vp["quality-q-vals"][node] = record.letter_annotations["phred_quality"]
    return graph


def set_edges(graph, read, threshold, stats):
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
                likelihood_calculation = likelihood_operations.get_alignment_likelihood(graph, read.cigartuples,
                                                               graph.vertex_properties["quality"][query_node],
                                                               stats)
                graph.ep["likelihood"][edge] = likelihood_calculation[0]
                stats = likelihood_calculation[1]
    return graph, stats


def save_and_draw_graph(graph, xml_out=None, figure_out=None, v_color=None, e_color=None, v_size=1):
    # default colors
    v_color_prop = "red"
    e_color_prop = "blue"

    if v_color:
        v_color_prop = graph.vertex_properties[v_color]
    if e_color:
        e_color_prop = graph.edge_properties[e_color]
    if xml_out:
        graph.save(xml_out)
    if figure_out:
        pos = sfdp_layout(graph)
        graph_draw(graph, vertex_fill_color=v_color_prop,
                   edge_color=e_color_prop,
                   pos=pos, vertex_size=v_size, output=figure_out)


def get_components(graph_obj, msg="", wildcards="", sub_dir=None, graph_xml=None, graph_figure=None, v_color=None, e_color=None):
    v_concom = graph_obj.vertex_properties["component-label"]
    all_components, hist = label_components(graph_obj, vprop=v_concom, directed=False)
    if msg:
        sys.stderr.write("\n\n{}:\n".format(msg))
        sys.stderr.write("Number of connected components: {}\n\n".format(max(all_components.a)))
        sys.stderr.write("Histogram of connected components:\n" + str(hist) + "\n\n")
        sys.stderr.write("connected components with more than one vertex are:\n\n")
    connected_components = []

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
                set_dir(sub_dir)
                xml = "{}/{}-{}.subgraph.xml.gz".format(sub_dir, wildcards, str(comp_nr))
                figure = "{}/{}-{}.subgraph.pdf".format(sub_dir, wildcards, str(comp_nr))
                save_and_draw_graph(conn_component, xml_out=xml, figure_out=figure, v_color=v_color, e_color=e_color, v_size=9)

    # optional output: all connected components of the graph:
    # set colors for nodes and edges and draw each component
    if graph_xml:
       save_and_draw_graph(graph_obj, xml_out=graph_xml)
    if graph_figure:
       save_and_draw_graph(graph_obj, figure_out=graph_figure, v_color=v_color, e_color=e_color)
    return connected_components




