import sys
import os
import gzip
import pysam
from Bio import SeqIO
from graph_tool.all import *

sys.stderr = open(snakemake.log[0], "w")

sam = snakemake.input.get("sam")
reads = snakemake.input.get("fastq")
threshold = snakemake.params.get("threshold", "")
mut_total = snakemake.params.get("mut_total", "")
mut_subst = snakemake.params.get("mut_subst", "")
n_snp = 0
mut_ins = snakemake.params.get("mut_ins", "")
n_ins = 0
mut_del = snakemake.params.get("mut_del", "")
n_del = 0
format = ""

if reads.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
    format = "fastq"
if reads.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")):
    format = "fasta"

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
graph.vertex_properties["seqence"] = v_seq
v_qual = graph.new_vertex_property("vector<float>")
graph.vertex_properties["quality"] = v_qual
e_dist = graph.new_edge_property("int")
graph.edge_properties["distance"] = e_dist
e_cs = graph.new_edge_property("string")  # cigar string from cs tag, short or long option can be selected in minimap2 rule
graph.edge_properties["cs-tag"] = e_cs
e_cig = graph.new_edge_property("string")  # cigar string from mandatory field 6 (sam format)
graph.edge_properties["cigar"] = e_cig
e_lh = graph.new_edge_property("float")
graph.edge_properties["likelihood"] = e_lh

node = graph.add_vertex()
v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=sample)
v_name[node] = ""
v_seq[node] = ""
v_qual[node] = ""
nodes = [node]
edges = []

# add vertices
def add_verticles(_reads, idx_nodes):
    for record in SeqIO.parse(_reads, format):
        idx_nodes += 1
        node = graph.add_vertex()
        v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=sample)
        v_name[node] = (record.id).split(' ', 1)[0]
        v_seq[node] = record.seq
        v_qual[node] = [10**(-1 * i / 10) for i in record.letter_annotations["phred_quality"]]
        nodes.append(node)

if reads.endswith(".gz"):
    with gzip.open(reads, "rt") as _reads:
        add_verticles(_reads, idx_nodes)
else:
    with open(reads, "rU") as _reads:
        add_verticles(_reads, idx_nodes)

# add edges from all-vs-all alignment of reads (please see rule minimap2)
sam = pysam.AlignmentFile(sam, "rb")
for read in sam.fetch(until_eof=True):
    if read.has_tag("NM") and not read.query_name == read.reference_name:
        nm = 0  # edit distance between query and reference sequence
        cig = ""  # elements of cigar string
        likelihood = 1.0
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
            if (not graph.vertex_index[node2] in graph.get_all_neighbors(node1)):
                edge = graph.add_edge(node1, node2)
                e_dist[edge] = nm
                e_cs[edge] = cig
                e_cig[edge] = read.cigarstring
                # pysam operations:
                    # M   BAM_CMATCH        0  -> 1 - basequal
                    # I   BAM_CINS          1  -> (1 - mutrate) * 1/3 * basequal + mutrate * (1 - basequal)
                    # D   BAM_CDEL          2  -> (1 - mutrate) * 1/3 * basequal + mutrate * (1 - basequal)
                    # N   BAM_CREF_SKIP     3  ??
                    # S   BAM_CSOFT_CLIP    4  ??
                    # H   BAM_CHARD_CLIP    5  ??
                    # P   BAM_CPAD          6  ??
                    # =   BAM_CEQUAL        7  -> 1 - basequal
                    # X   BAM_CDIFF         8  -> (1 - mutrate) * 1/3 * basequal + mutrate * (1 - basequal)
                    # B   BAM_CBACK         9  ??
                qual_idx = 0
                for (op, length) in read.cigartuples:
                    mutrate = mut_total
                    n1 = graph.vertex_properties["quality"][node1]
                    n2 = graph.vertex_properties["quality"][node2]
                    if op == 7 or op == 0:  # on match
                        for i in n1[qual_idx:length-1]:
                            likelihood *= 1 - i
                        for i in n2[qual_idx:length - 1]:
                            likelihood *= 1 - i
                        qual_idx = length
                    if op == 8 or op == 1 or op == 2:
                        if op == 8 and mut_subst:  # on mismatch: substitution/snp
                            mutrate = mut_subst
                            n_snp += length
                        if op == 1 and mut_ins:  # on mismatch: insertion
                            mutrate = mut_ins
                            n_ins += length
                        if op == 2 and mut_del:  # on mismatch: deletion
                            mutrate = mut_del
                            n_del += length
                        for i in n1[qual_idx:length-1]:
                            likelihood *= (1 - mutrate) * float(1/3) * i + mutrate * (1 - i)  # or (1-basequal) * 1/3 * basequal + mutrate * (1-basequal) ?
                        for i in n2[qual_idx:length - 1]:
                            likelihood *= (1 - mutrate) * float(1/3) * i + mutrate * (1 - i)  # or (1-basequal) * 1/3 * basequal + mutrate * (1-basequal) ?
                        qual_idx = length
                e_lh[edge] = likelihood
                edges.append(edge)
sam.close()

n_nodes = len(nodes)
n_edges = len(edges)
sys.stderr.write("graph construction summary for sample {}:\n nodes:\t{}\n edges:\t{}\n max distance:\t{}".format(sample, n_nodes, n_edges, max_NM))
sys.stderr.write("\n\tmutations:\n\t\tsnp's and substitutions: {}\n\t\tinsertions: {}\n\t\tdeletions: {}\n\n".format(n_snp, n_ins, n_del))

graph.save(snakemake.output.get("graph_xml"))
pos = sfdp_layout(graph)
graph_draw(graph, vertex_color=[1, 1, 1, 0],
           edge_color=graph.edge_properties["likelihood"],
           pos=pos, vertex_size=1, output=snakemake.output.get("graph_figure"))
# hist = distance_histogram(graph, weight=e_lh, samples=n_nodes)
# print(hist)
