import sys
import os
from graph_tool.all import *
import pulp

sys.stderr = open(snakemake.log[0], "w")

graph = load_graph(snakemake.input[0])
e_lh = graph.edge_properties["likelihood"]
v_id = graph.vertex_properties["id"]
v_concom = graph.new_vertex_property("int32_t")
graph.vertex_properties["component-label"] = v_concom

# extract connected components
all_components, hist = label_components(graph, vprop=v_concom)
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
        # ### for testing and debugging:
        # pos = sfdp_layout(conn_component)
        # sub_dir = "/home/tarja/Schreibtisch/workspace_bioinformatik/NodeRAD/.test/results/subgraphes/{}".format(snakemake.wildcards.get('sample'))
        # if not os.path.exists(sub_dir):
        #     os.makedirs(sub_dir)
        # graph_draw(conn_component, vertex_fill_color=conn_component.vertex_properties['component-label'],
        #            edge_color=conn_component.edge_properties["likelihood"],
        #            pos=pos, vertex_size=9, output=sub_dir+"/"+str(comp_nr)+"-view_subgraph.pdf")
        # conn_component.save(sub_dir+"/"+str(comp_nr)+"-subgraph.xml.gz")
    else:
        sporadic_vertices.append(conn_component)

# ### for testing and debugging:
# sub_dir = "/home/tarja/Schreibtisch/workspace_bioinformatik/NodeRAD/.test/results/subgraphes/{}".format(snakemake.wildcards.get('sample'))
# if not os.path.exists(sub_dir):
#     os.makedirs(sub_dir)
# graph.save(sub_dir+"/"+"component.xml.gz")
# pos = sfdp_layout(graph)
# graph_draw(graph, vertex_fill_color=graph.vertex_properties['component-label'],
#            # edge_color=graph.edge_properties["likelihood"],
#            pos=pos, vertex_size=1, output=sub_dir+"/"+"components.pdf")

# ILP on each connected component


# # function:
# # representatives_ilp = pulp.LpProblem("Find representatives for loci", pulp.LpMaximize)
#
# # restrictions:


for comp in connected_components:
    for node in comp.vertices():
        print(comp.vertex_properties['id'][node])
    print("--------------------------------------------------")

######## DEBUG ########
with open(snakemake.output[0], 'w') as f:
    print("success!", file=f)
#######################
