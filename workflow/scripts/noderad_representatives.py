import sys
from graph_tool.all import *
import pulp
import math
import noderad_solver_wrapper
import noderad_graph_operations


sys.stderr = open(snakemake.log[0], "w")


### input
graph = load_graph(snakemake.input[0])


### required output
repres = snakemake.output.get("representatives", "")
repres_xml = snakemake.output.get("representatives_xml", "")
repres_figure = snakemake.output.get("repesentatives_figure", "")


### optional output
connected_components_xml = snakemake.output.get("connected_components_xml", "")
connected_components_figure = snakemake.output.get("connected_components_figure", "")
dir_subgraphs = snakemake.output.get("dir_subgraphs", "")


### solver configuration wrapper for PuLP:
solver = snakemake.params.get("solver", None)
mip = snakemake.params.get("mip", None)
timelimit = snakemake.params.get("timelimit", None)
gaprel = snakemake.params.get("gaprel", None)
gapabs = snakemake.params.get("gapabs", None)
maxnodes = snakemake.params.get("maxnodes", None)
maxmemory = snakemake.params.get("maxmemory", None)
threads = snakemake.threads
spec_solver = ""

# apply configurable solver options, if not specified, the default solver is used
if solver:
    if solver in pulp.apis.list_solvers(onlyAvailable=True):
        spec_solver = noderad_solver_wrapper.check_params(solver, mip, timelimit, gaprel, gapabs, maxnodes, maxmemory, threads)
    else:
        sys.stderr.write(
            'The selected solver is not supported or is not available. The calculation is continued with the default solvers of Pulp (CBC and CHOCO).')


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
    likelihood_log_sum = 0
    likelihood_sum = 0
    for neighbor in graph.get_in_neighbors(graph.vertex_index[node]):
        likelihood = graph.edge_properties['likelihood'][graph.edge(neighbor, node)]
        likelihood_log_sum += math.log10(likelihood)
        likelihood_sum += likelihood
    v_log_sum_in[node] = likelihood_log_sum
    v_sum_in[node] = likelihood_sum


### extract connected components
message = "CONNECTED COMPONENTS based on the graph construction from the edit distances (minimap2)"
connected_components = noderad_graph_operations.get_components(graph, message, snakemake.wildcards.get('sample'),
                                                                   dir_subgraphs, connected_components_xml,
                                                                   connected_components_figure, v_color="component-label",
                                                                   e_color="likelihood")


### ILP on each connected component

with open(repres, 'a+') as f:
    print("{}\t{}\t{}\t{}".format("Component-ID", "Total number of representatives", "ILP-solution", "Node-ID of representatives"), file=f)

for (comp, comp_nr) in connected_components:
    sys.stderr.write("\n\nSolutions for each number of representatives for connected component {}:\n".format(comp_nr))

    # data:
    nodes = list([comp.vertex_index[node] for node in comp.vertices()])

    # determining the lower limit for the solvability of the ILP
    comp_size = len(nodes)
    vertex_degrees = []
    for node in nodes:
        vertex_degrees.append(comp.vertex(node).in_degree() + comp.vertex(node).out_degree())
    vertex_degrees.sort(reverse=True)
    minimum_representatives = noderad_graph_operations.get_minimum_representatives(vertex_degrees, comp_size)

    for n in range(minimum_representatives, comp_size + 1):
        # decision variables
        r = pulp.LpVariable.dicts("representatives", nodes, 0, 1, pulp.LpBinary)
        z = pulp.LpVariable.dicts("indicator", (nodes, nodes), 0, 1, pulp.LpBinary)

        # function:
        model_representatives = pulp.LpProblem("Find_representatives_for_loci", pulp.LpMaximize)
        model_representatives += pulp.lpSum([comp.vertex_properties["log-sum-in-edges"][node] * r[node] for node in nodes])

        # restrictions:
        for node in nodes:
            for neighbor in comp.get_in_neighbors(comp.vertex_index[node]):
                # z[i][j] = 1 <=> i is adjacent to j AND r[j] == 1
                model_representatives += (z[neighbor][node] == r[node]), "{i}_adjacent_to_{j}_and_r[{j}]_is_1".format(
                    i=node, j=neighbor)

            # each node has a representative or is a representative itself
            model_representatives += pulp.lpSum(
                [z[node][neighbor] for neighbor in comp.get_out_neighbors(comp.vertex_index[node])]) >= (1 - r[
                node]), "each_node[{i}]_is_representative_or_has_representative[{j}]".format(i=node, j=neighbor)

        # n representatives must be found
        model_representatives += pulp.lpSum([r[node] for node in nodes]) == n, "for_each_number_of_representatives"

        if spec_solver:
            model_representatives.solve(spec_solver)
        else:
            model_representatives.solve()

        # write results to graph properties
        n_representants = 0
        if pulp.LpStatus[model_representatives.status] == 'Optimal':
            for var in model_representatives.variables():
                if (var.name).startswith('indicator'):
                    if var.varValue == 1:
                        idx_1 = var.name.split("_")[-2]
                        idx_2 = var.name.split("_")[-1]
                        node_1 = find_vertex(graph, v_id, comp.vertex_properties["id"][idx_1])[0]
                        node_2 = find_vertex(graph, v_id, comp.vertex_properties["id"][idx_2])[0]
                        edg = graph.edge(node_1, node_2)
                        e_n_indicator[edg].append(n)

                if (var.name).startswith('representatives'):
                    if var.varValue == 1:
                        n_representants += 1
                        idx = var.name.split("_")[-1]
                        node = find_vertex(graph, v_id, comp.vertex_properties["id"][idx])[0]
                        solution = pulp.value(model_representatives.objective.value())
                        v_n_repres[node].append(n)
                        v_sol_repres[node].append(solution)

                        with open(repres, 'a+') as f:
                            print("{}\t{}\t{}\t{}".format(comp_nr, n, solution, graph.vertex_properties["id"][node]), file=f)

        sys.stderr.write("\n\tstatus: {}".format(pulp.LpStatus[model_representatives.status]))
        sys.stderr.write("\n\tFor n = {} were {} representatives "
                         "found with solution value: {}.\n".format(n, n_representants,
                                                                   pulp.value(model_representatives.objective.value())))

### save graph with new properties
graph.save(repres_xml)

### draw graph
# colors for nodes: frequency of being selected as representative,
# colors for edges: frequency of being selected as belonging to a representative
v_color = graph.new_vertex_property('int')
graph.vertex_properties['repr_color'] = v_color
for node in graph.vertices():
    repr_n = graph.vertex_properties['total-number-of-representatives'][node]
    if repr_n:
        v_color[node] = len(repr_n)
    else:
        v_color[node] = 0
e_color = graph.new_edge_property('int')
graph.edge_properties['indicator_color'] = e_color
for edge in graph.edges():
    indicator_n = graph.edge_properties['indicator'][edge]
    if indicator_n:
        e_color[edge] = len(indicator_n)
    else:
        e_color[edge] = 0

if repres_figure:
    pos = sfdp_layout(graph)
    graph_draw(graph, vertex_fill_color=graph.vertex_properties['repr_color'],
               edge_color=graph.edge_properties['indicator_color'],
               pos=pos, vertex_size=1, output=repres_figure)
