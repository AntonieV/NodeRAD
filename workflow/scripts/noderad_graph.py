import sys
import os
import gzip
import pysam
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from graph_tool.all import *

sys.stderr = open(snakemake.log[0], "w")

sam = snakemake.input.get("sam")
reads = snakemake.input.get("fastq")
threshold = snakemake.params.get("threshold", "")

# default value of threshold for maximum distance value
if threshold == "":
    threshold = 23
else:
    threshold = int(threshold)

sample = os.path.splitext(os.path.basename(reads))[0]  # name of sample/reads
if reads.endswith(".gz"):
    sample = os.path.splitext(sample)[0]

# init graph
idx_nodes = 0  # index of nodes
max_NM = 0     # highest edit distance
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
e_cig = graph.new_edge_property("string")
graph.edge_properties["cig"] = e_cig

node = graph.add_vertex()
v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=sample)
v_name[node] = ""
v_seq[node] = ""
v_qual[node] = ""
nodes = [node]
edges = []

# add vertices
def add_verticles(_reads, idx_nodes):
    for (id, seq, qual) in FastqGeneralIterator(_reads):
        idx_nodes += 1
        node = graph.add_vertex()
        v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=sample)
        v_name[node] = (id).split(' ', 1)[0]
        v_seq[node] = seq
        v_qual[node] = qual
        nodes.append(node)

if reads.endswith(".gz"):
    with gzip.open(reads, "rt") as _reads:
        add_verticles(_reads, idx_nodes)
else:
    with open(reads, "rU") as _reads:
        add_verticles(_reads, idx_nodes)

# add edges
sam = pysam.AlignmentFile(sam, "rb")
for read in sam.fetch(until_eof=True):
    if read.has_tag("NM") and not read.query_name == read.reference_name:
        nm = 0
        cig = ""
        for (tag, tag_val) in read.get_tags():
            if tag == "NM":
                nm = tag_val
            if tag == "cs":
                cig = tag_val
        if nm <= threshold:
            if max_NM < nm:
                max_NM = nm
            node1 = find_vertex(graph, v_name, read.query_name)[0]
            node2 = find_vertex(graph, v_name, read.reference_name)[0]
            edge = graph.add_edge(node1, node2)
            e_dist[edge] = nm
            e_cig[edge] = cig
            edges.append(edges)
sam.close()

n_nodes = len(nodes)
n_edges = len(edges)
sys.stderr.write("graph construction summary for sample {}:\n nodes:\t{}\n edges:\t{}\n max distance:\t{}".format(sample, n_nodes, n_edges, max_NM))

graph.save(snakemake.output.get("graph_xml"))
pos = sfdp_layout(graph)
graph_draw(graph, vertex_color=[1, 1, 1, 0],
           edge_color=graph.edge_properties["dist"], pos=pos, vertex_size=1,
           output=snakemake.output.get("graph_figure"))

