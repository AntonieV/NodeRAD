from graph_tool.all import *
import sys
import os

def get_properties(graphtype):
    if graphtype == "distance-graph":
        return [(name, prop) for (name, prop, prop_type) in get_new_properties("init-graph")]

    # if graphtype == "optimized-graph":  #ToDo


def get_new_properties(graphtype):
    if graphtype == "init-graph":
        return [("v_id", "id", "string"), ("v_name", "name", "string"),
                ("v_seq", "sequence", "string"), ("v_qual", "quality", "vector<float>"),
                ("e_dist", "distance", "int"),
                ("e_cs", "cs-tag", "string"),  # cigar string from cs tag, short or long option can be selected in minimap2 rule
                ("e_cig", "cigar", "string"),  # cigar string from mandatory field 6 (sam format)
                ("e_lh", "likelihood", "float")]

    if graphtype == "distance-graph":
        ### new properties for ILP calculation
        # v_concom: index number of the connected component to which the node belongs to
        # v_sum_in: sum of likelihood of all incoming edges of the node
        # v_n_repres: list of all n values where the node was chosen as representative for the optimal solution at the ilp
        # v_sol_repres: list of all ilp solutions where the node is representative, the corresponding n is located in
        # v_n_repres at the same index. Therefore (v_n_repres[i], v_sol_repres[i]) defines the optimal solution
        # if v_n_repres[i] representants are used
        # e_n_indicator: list of total number of representatives where the edge (source, target) is set on z[source][target] == 1
        return [("v_concom", "component-label", "int32_t"), ("v_sum_in", "sum-in-edges", "double"),
                ("v_n_repres", "total-number-of-representatives", "vector<long>"),
                ("v_sol_repres", "ilp-solution", "vector<double>"), ("e_n_indicator", "indicator", "vector<long>")]

    # if graphtype == "optimized-graph":  #ToDo


def get_components(graph_obj, msg="", wildcards="", sub_dir=None, graph_xml=None, graph_figure=None):
    v_concom = graph_obj.vertex_properties["component-label"]
    all_components, hist = label_components(graph_obj, vprop=v_concom)
    sys.stderr.write("\n\n{}:\n".format(msg))
    sys.stderr.write("Number of connected components: {}\n\n".format(max(all_components.a)))
    sys.stderr.write("Histogram of connected components:\n" + str(hist) + "\n\n")
    sys.stderr.write("connected components with more than one vertex are:\n\n")
    connected_components = []
    sporadic_vertices = []

    for comp_nr in range(max(all_components)):
        # a is a shortcut for the get_array() method as an attribute,
        # please see https://graph-tool.skewed.de/static/doc/graph_tool.html#graph_tool.PropertyMap.get_array
        conn_component = GraphView(graph_obj, vfilt=all_components.a == comp_nr)
        conn_component = Graph(conn_component, prune=True)
        if hist[comp_nr] > 1:
            connected_components.append(tuple([conn_component, comp_nr]))
            sys.stderr.write("\t" + str(conn_component) + "\n")

            # optional output files for each subgraph:
            if sub_dir:
                if not os.path.exists(sub_dir):
                    os.makedirs(sub_dir)
                conn_component.save(
                    "{}/{}-{}.subgraph.xml.gz".format(sub_dir, wildcards, str(comp_nr)))
                pos = sfdp_layout(conn_component)
                graph_draw(conn_component, vertex_fill_color=conn_component.vertex_properties['component-label'],
                           edge_color=conn_component.edge_properties["likelihood"],
                           pos=pos, vertex_size=9,
                           output="{}/{}-{}.subgraph.pdf".format(sub_dir, wildcards,
                                                                 str(comp_nr)))
        else:
            sporadic_vertices.append(conn_component)

    # optional output: connected components of the graph:
    if graph_xml:
        graph_obj.save(graph_xml)
    if graph_figure:
        pos = sfdp_layout(graph_obj)
        graph_draw(graph_obj, vertex_fill_color=graph_obj.vertex_properties['component-label'],
                   edge_color=graph_obj.edge_properties["likelihood"],
                   pos=pos, vertex_size=1, output=graph_figure)
    return connected_components


def get_minimum_representatives(vertex_degrees, comp_size):
    minimum_degrees = 0
    required_nodes = 0
    for i in vertex_degrees:
        minimum_degrees += i
        required_nodes += 1
        if minimum_degrees >= comp_size:
            return required_nodes
    sys.exit("The tested component is not a connected component.")
