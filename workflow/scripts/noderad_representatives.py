import sys
import os
from graph_tool.all import *
import pulp
import math
import numpy as np

sys.stderr = open(snakemake.log[0], "w")

graph = load_graph(snakemake.input[0])
e_lh = graph.edge_properties["likelihood"]
v_id = graph.vertex_properties["id"]
v_concom = graph.new_vertex_property("int32_t")
graph.vertex_properties["component-label"] = v_concom

# optional output:
connected_components_xml = snakemake.output.get("connected_components_xml", "")
connected_components_figure = snakemake.output.get("connected_components_figure", "")
dir_subgraphs = snakemake.output.get("dir_subgraphs", "")

# extract connected components
all_components, hist = label_components(graph, vprop=v_concom)
# a is a shortcut for the get_array() method as an attribute, please see https://graph-tool.skewed.de/static/doc/graph_tool.html#graph_tool.PropertyMap.get_array
sys.stderr.write("number of connected components: {}\n\n".format(max(all_components.a)))
sys.stderr.write("histogram of connected components:\n"+str(hist)+"\n\n")
connected_components = []
sporadic_vertices = []
sys.stderr.write("connected components with more than one vertex are:\n\n")
for comp_nr in range(max(all_components)):
    conn_component = GraphView(graph, vfilt=all_components.a == comp_nr)
    conn_component = Graph(conn_component, prune=True)
    if hist[comp_nr] > 1:
        connected_components.append(conn_component)
        sys.stderr.write("\t" + str(conn_component) + "\n")

        # optional output files for each subgraph:
        if dir_subgraphs:
            if not os.path.exists(dir_subgraphs):
                os.makedirs(dir_subgraphs)
            conn_component.save("{}/{}-{}.subgraph.xml.gz".format(dir_subgraphs, snakemake.wildcards.get('sample'), str(comp_nr)))
            pos = sfdp_layout(conn_component)
            graph_draw(conn_component, vertex_fill_color=conn_component.vertex_properties['component-label'],
                       edge_color=conn_component.edge_properties["likelihood"],
                       pos=pos, vertex_size=9, output="{}/{}-{}.subgraph.pdf".format(dir_subgraphs, snakemake.wildcards.get('sample'), str(comp_nr)))

    else:
        sporadic_vertices.append(conn_component)

# optional: connected components of the graph:
if connected_components_xml:
    graph.save(connected_components_xml)
if connected_components_figure:
    pos = sfdp_layout(graph)
    graph_draw(graph, vertex_fill_color=graph.vertex_properties['component-label'],
               edge_color=graph.edge_properties["likelihood"],
               pos=pos, vertex_size=1, output=connected_components_figure)

# ILP on each connected component

for comp in connected_components:
    # data:
    nodes = list([comp.vertex_index[node] for node in comp.vertices()])
    print(nodes)
    edges = list([(comp.vertex_index[ed.source()], comp.vertex_index[ed.target()]) for ed in comp.edges()])

    for n in range(1, len(nodes) + 1):
        # decision variables
        r = pulp.LpVariable.dicts("representatives", nodes, 0, 1, pulp.LpBinary)
        z = pulp.LpVariable.dicts("indicator", (nodes, nodes), 0, 1, pulp.LpBinary)

        # function:
        model_representatives = pulp.LpProblem("Find_representatives_for_loci", pulp.LpMaximize)
        # includes the constrait z[i][j] = 1 :
            # get_all_edges(node) guarantees: z[i][j] = 1 <=> i is adjacent to j AND
            # r == 1 guarantees: z[i][j] = 0 if r[i] == 0
        model_representatives += pulp.lpSum([z[source][target] * math.log(comp.edge_properties['likelihood'][comp.edge(source, target)]) for node in nodes for (source, target) in comp.get_all_edges(comp.vertex_index[node]) if r[source] == 1])

        # # restrictions:
        model_representatives += pulp.lpSum([z[source][target] for node in nodes for (source, target) in comp.get_all_edges(comp.vertex_index[node])]) >= 1, "each_node_is_representative_or_has_representative"
        model_representatives += pulp.lpSum([r[node] for node in nodes]) == n, "for_each_number_of_representatives"
        model_representatives.solve()
        print("Best solution for {} representatives is {} with: ".format(n, pulp.value(model_representatives.objective.value())))
        for var in model_representatives.variables():
            if var.varValue == 1:
                print("{} = {}".format(var.name, var.varValue))
        with open(snakemake.output.get("representatives"), 'a+') as f:
            print("status: ", pulp.LpStatus[model_representatives.status], file=f)
            print("objective-value: ", pulp.value(model_representatives.objective.value()), file=f)
        # for var in model_representatives.variables():
        #     with open(snakemake.output.get("representatives"), 'a+') as f:
        #         if var.varValue == 1:
        #             print("Best solution for {} representatives is {}".format())
        #         print("{} = {}".format(var.name, var.varValue), file=f)
