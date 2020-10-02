import sys
import os
import mappy as mp  # Minimap2 library (Python binding), see https://github.com/lh3/minimap2/tree/master/python
from graph_tool.all import *

sys.stderr = open(snakemake.log[0], "w")

input = snakemake.input[0]
run_name = os.path.splitext(os.path.basename(input))[0]  # name of sample/input
if input.endswith(".gz"):
    run_name = os.path.splitext(run_name)[0]

alignment = mp.Aligner(input)  # load or build index
if not alignment:
    raise Exception("ERROR: Failed to build alignment index.")

thread_buffer = mp.ThreadBuffer()
reads = mp.fastx_read(input)

# init graph
idx_nodes = 0  # index of nodes
n_nodes = len(alignment.seq_names)  # number of nodes
n_edges = 0  # number of edges
max_NM = 0  # highest edit distance
graph = Graph(directed=False)
v_id = graph.new_vertex_property("string")
graph.vertex_properties["id"] = v_id
v_name = graph.new_vertex_property("string")
graph.vertex_properties["name"] = v_name
v_seq = graph.new_vertex_property("string")
graph.vertex_properties["seq"] = v_seq
v_qual = graph.new_vertex_property("string")
graph.vertex_properties["qual"] = v_qual
e_dist = graph.new_edge_property("int")
graph.edge_properties["dist"] = e_dist

node = graph.add_vertex()
v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=run_name)
v_name[node] = ""
v_seq[node] = ""
v_qual[node] = ""
nodes = [node]

# add vertices

for name, seq, qual in mp.fastx_read(snakemake.input[0]):
    idx_nodes += 1
    node = graph.add_vertex()
    v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=run_name)
    v_name[node] = name
    v_seq[node] = seq
    v_qual[node] = qual
    nodes.append(node)
    # print(graph.get_vertices(node))

# for name, seq, qual in mp.fastx_read(snakemake.input[0]):  # read a fasta/q sequence

for node in nodes:
    for hit in alignment.map(v_seq[node], buf=thread_buffer):  # traverse alignments
        #     for hit in alignment.map(seq, buf=thread_buffer):
        if not v_name[node] == hit.ctg:
            # if not name == hit.ctg:
            node2 = find_vertex(graph, v_name, hit.ctg)
            if len(node2) == 1:
                node2 = node2[0]
                print(node2)
            else:
                sys.exit(
                    "There are ambiguous read signatures."
                )
            graph.add_edge(node, node2)
            n_edges += 1
            # print("\n\n" + str(hit) + "\nName: " + name + " Seq:  " + seq + " Qual: " + qual)  ###
            # print("{}\t{}\t{}\t{}\t{}".format(hit.mlen, hit.ctg, hit.NM, hit.cs, hit.cigar_str)) ###
            # print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(hit.mlen, hit.ctg, hit.r_st, hit.r_en, hit.cigar_str, hit.cs, hit.NM)) ###
            # alignment = next(aligner.map(seq)) ###
            if max_NM < hit.NM:
                max_NM = hit.NM

print("samples: {}; hits: {}; max: {}".format(n_nodes, n_edges, max_NM))  ###

############################### Debugging - remove later ###############################
# s = a.seq("MT_human", 100, 200)     ### retrieve a subsequence from the index
# print(mp.revcomp(s))                ### reverse complement
x = list(mp.fastx_read(snakemake.input[0]))[3][0]  ### indexing on type: generator ###
print(x)  ###
print(v_id[nodes[1]])
# print(graph.list_properties(graph.get_vertices()["name"]))
with open(snakemake.output[0], "w") as out:  ###
    print("samples: {}; hits: {}; max: {}".format(n_nodes, n_edges, max_NM), file=out)  ###

########################################################################################
graph_out_dir = "results/graph"
if not os.path.exists(graph_out_dir):
    os.mkdir(graph_out_dir)

graph.save("{out}/{sample}.xml.gz".format(out=graph_out_dir, sample=run_name))

graph_draw(graph, output="{out}/{sample}.pdf".format(out=graph_out_dir, sample=run_name))
