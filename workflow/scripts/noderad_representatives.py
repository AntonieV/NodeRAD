import sys
from graph_tool.all import *
import pulp
import math
import noderad_solver_wrapper
import noderad_connected_components

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


### vertex properties
v_id = graph.vertex_properties["id"]
v_name = graph.vertex_properties["name"]
v_seq = graph.vertex_properties["sequence"]
v_qual = graph.vertex_properties["quality"]

### edge properties
e_dist = graph.edge_properties["distance"]
e_cs = graph.edge_properties["cs-tag"]
e_cig = graph.edge_properties["cigar"]
e_lh = graph.edge_properties["likelihood"]


### additional vertex and edge properties for graph
# index number of the connected component to which the node belongs to
v_concom = graph.new_vertex_property("int32_t")
graph.vertex_properties["component-label"] = v_concom

# sum of likelihood of all incoming edges of the node
v_sum_in = graph.new_vertex_property("double")
graph.vertex_properties["sum-in-edges"] = v_sum_in

# array of all n values where the node was chosen as representative for the optimal solution at the ilp
v_n_repres = graph.new_vertex_property("vector<long>")
graph.vertex_properties["total-number-of-representatives"] = v_n_repres

# array of all ilp solutions where the node is representative, the corresponding n is located in
# v_n_repres at the same index
# therefore (v_n_repres[i], v_sol_repres[i]) defines the optimal solution if v_n_repres[i] representants are used
v_sol_repres = graph.new_vertex_property("vector<double>")
graph.vertex_properties["ilp-solution"] = v_sol_repres

# array of total number of representatives where the edge (source, target) is set on z[source][target] == 1
e_n_indicator = graph.new_edge_property("vector<long>")
graph.edge_properties["indicator"] = e_n_indicator


### pre-calculation of the sum of likelihoods of incoming edges for each node (to improve the runtime of the ILP)
for node in graph.vertices():
    likelihood_sum = 0
    for neighbor in graph.get_in_neighbors(graph.vertex_index[node]):
        likelihood_sum += math.log10(graph.edge_properties['likelihood'][graph.edge(neighbor, node)])
    v_sum_in[node] = likelihood_sum


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


### extract connected components
message = "CONNECTED COMPONENTS based on the graph construction from the edit distances (minimap2)"
connected_components = noderad_connected_components.get_components(graph, message, snakemake.wildcards.get('sample'),
                                                                   dir_subgraphs, connected_components_xml,
                                                                   connected_components_figure)


### ILP on each connected component

with open(repres, 'a+') as f:
    print("{}\t{}\t{}".format("Component-ID", "Total number of representatives", "Node-ID of representatives"), file=f)

for (comp, comp_nr) in connected_components:
    sys.stderr.write("\n\nSolutions for each number of representatives for connected component {}:\n".format(comp_nr))
    # data:
    nodes = list([comp.vertex_index[node] for node in comp.vertices()])

    for n in range(1, len(nodes) + 1):
        # decision variables
        r = pulp.LpVariable.dicts("representatives", nodes, 0, 1, pulp.LpBinary)
        z = pulp.LpVariable.dicts("indicator", (nodes, nodes), 0, 1, pulp.LpBinary)

        # function:
        model_representatives = pulp.LpProblem("Find_representatives_for_loci", pulp.LpMaximize)
        # model_representatives += pulp.lpSum([z[neighbor][node] * math.log10(comp.edge_properties['likelihood'][comp.edge(neighbor, node)]) for node in nodes for neighbor in comp.get_in_neighbors(comp.vertex_index[node])])
        model_representatives += pulp.lpSum([comp.vertex_properties["sum-in-edges"][node] * r[node] for node in nodes])

        # restrictions:
        for node in nodes:
            for neighbor in comp.get_in_neighbors(comp.vertex_index[node]):
                # z[i][j] = 1 <=> i is adjacent to j AND r[j] == 1
                model_representatives += (z[neighbor][node] == r[node]), "{i}adjacent_to_{j}_and_r[{j}]_is_1".format(
                    i=node, j=neighbor)

            # each node has a representative or is a representative itself
            model_representatives += pulp.lpSum(
                [z[node][neighbor] for neighbor in comp.get_out_neighbors(comp.vertex_index[node])]) >= (1 - r[
                node]), "each_node[{i}]_is_representative_or_has_representative[{j}]".format(i=node, j=neighbor)

        model_representatives += pulp.lpSum([r[node] for node in nodes]) == n, "for_each_number_of_representatives"
        # model_representatives.writeLP(snakemake.output.get("representatives"))

        if spec_solver:
            model_representatives.solve(spec_solver)
        else:
            model_representatives.solve()

        # write results to graph properties
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
                        idx = var.name.split("_")[-1]
                        node = find_vertex(graph, v_id, comp.vertex_properties["id"][idx])[0]
                        v_n_repres[node].append(n)
                        v_sol_repres[node].append(pulp.value(model_representatives.objective.value()))

                        with open(repres, 'a+') as f:
                            print("{}\t{}\t{}".format(comp_nr, n, graph.vertex_properties["id"][node]), file=f)

        # with open(snakemake.output.get("representatives"), 'a+') as f:
        #     print("status: ", pulp.LpStatus[model_representatives.status], file=f)
        #     print("Best solution for {} representatives is {} with: ".format(n, pulp.value(
        #         model_representatives.objective.value())), file=f)
        sys.stderr.write("\n\tstatus: {}".format(pulp.LpStatus[model_representatives.status]))
        sys.stderr.write("\n\tBest solution for {} representatives is {}.\n".format(n, pulp.value(
                model_representatives.objective.value())))

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

pos = sfdp_layout(graph)
graph_draw(graph, vertex_fill_color=graph.vertex_properties['repr_color'],
           edge_color=graph.edge_properties['indicator_color'],
           pos=pos, vertex_size=1, output=repres_figure)

# for edge in edges, if len(e_n_indicator)==0 -> remove edge
