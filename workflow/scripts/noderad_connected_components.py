from graph_tool.all import *
import sys
import os

def get_components(graph_obj, msg, wildcards, sub_dir, graph_xml, graph_figure):
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

    # optional: connected components of the graph:
    if graph_xml:
        graph_obj.save(graph_xml)
    if graph_figure:
        pos = sfdp_layout(graph_obj)
        graph_draw(graph_obj, vertex_fill_color=graph_obj.vertex_properties['component-label'],
                   edge_color=graph_obj.edge_properties["likelihood"],
                   pos=pos, vertex_size=1, output=graph_figure)
    return connected_components
