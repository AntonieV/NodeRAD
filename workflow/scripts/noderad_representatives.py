import sys
import os
from graph_tool.all import *
import pulp
import math
import collections
# import pandas as pd
import solver_wrapper

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

# solver configuration wrapper for PuLP:
solver = snakemake.params.get("solver", None)
mip = snakemake.params.get("mip", None)
timelimit = snakemake.params.get("timelimit", None)
gaprel = snakemake.params.get("gaprel", None)
gapabs = snakemake.params.get("gapabs", None)
maxnodes = snakemake.params.get("maxnodes", None)
maxmemory = snakemake.params.get("maxmemory", None)
threads = snakemake.threads
spec_solver = ""

if solver:
    if solver in pulp.apis.list_solvers(onlyAvailable=True):
        spec_solver = solver_wrapper.check_params(solver, mip, timelimit, gaprel, gapabs, maxnodes, maxmemory, threads)
    else:
        sys.stderr.write('The selected solver is not supported or is not available. The calculation is continued with the default solvers of Pulp (CBC and CHOCO).')

# extract connected components
all_components, hist = label_components(graph, vprop=v_concom)
# a is a shortcut for the get_array() method as an attribute, please see https://graph-tool.skewed.de/static/doc/graph_tool.html#graph_tool.PropertyMap.get_array
sys.stderr.write("Number of connected components: {}\n\n".format(max(all_components.a)))
sys.stderr.write("Histogram of connected components:\n"+str(hist)+"\n\n")
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

    for n in range(1, len(nodes) + 1):
        # decision variables
        r = pulp.LpVariable.dicts("representatives", nodes, 0, 1, pulp.LpBinary)
        z = pulp.LpVariable.dicts("indicator", (nodes, nodes), 0, 1, pulp.LpBinary)

        # function:
        model_representatives = pulp.LpProblem("Find_representatives_for_loci", pulp.LpMaximize)
        model_representatives += pulp.lpSum([z[neighbor][node] * math.log10(comp.edge_properties['likelihood'][comp.edge(neighbor, node)]) for node in nodes for neighbor in comp.get_in_neighbors(comp.vertex_index[node])])

        # restrictions:
        for node in nodes:
            for neighbor in comp.get_in_neighbors(comp.vertex_index[node]):
                # z[i][j] = 1 <=> i is adjacent to j AND r[j] == 1
                model_representatives += (z[neighbor][node] == r[node]), "{i}adjacent_to_{j}_and_r[{j}]_is_1".format(i=node, j=neighbor)

        for node in nodes:
            model_representatives += pulp.lpSum([z[node][neighbor] for neighbor in comp.get_out_neighbors(comp.vertex_index[node])]) >= 1, "each_node[{i}]_is_representative_or_has_representative[{j}]".format(i=node, j=neighbor)

        model_representatives += pulp.lpSum([r[node] for node in nodes]) == n, "for_each_number_of_representatives"
        # model_representatives.writeLP(snakemake.output.get("representatives"))

        if spec_solver:
            model_representatives.solve(spec_solver)
        else:
            model_representatives.solve()

        with open(snakemake.output.get("representatives"), 'a+') as f:
            print("status: ", pulp.LpStatus[model_representatives.status], file=f)
            print("Best solution for {} representatives is {} with: ".format(n, pulp.value(model_representatives.objective.value())), file=f)
