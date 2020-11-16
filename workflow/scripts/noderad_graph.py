import sys
import os
import gzip
import pysam
from Bio import SeqIO
from graph_tool.all import *
import noderad_graph_operations

sys.stderr = open(snakemake.log[0], "w")

### input
sam = snakemake.input.get("sam")
reads = snakemake.input.get("fastq")

### optional output
output_figure = snakemake.output.get("graph_figure", "")

### params for mutation rates and threshold
threshold = snakemake.params.get("threshold", "")
mut_total = snakemake.params.get("mut_total", "")
mut_subst = snakemake.params.get("mut_subst", "")
mut_ins = snakemake.params.get("mut_ins", "")
mut_del = snakemake.params.get("mut_del", "")

n_snp = 0      # counter for total number of substitutions/snp's
n_ins = 0      # counter for total number of insertions
n_del = 0      # counter for total number of deletions
idx_nodes = 0  # index of nodes
max_NM = 0     # highest edit distance
format = ""    # format of reads


### get reads format
if reads.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
    format = "fastq"
if reads.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")):
    format = "fasta"

### default value of threshold for maximum distance value
if threshold == "":
    threshold = 23
else:
    threshold = int(threshold)

sample = os.path.splitext(os.path.basename(reads))[0]  # name of sample/reads
if reads.endswith(".gz"):
    sample = os.path.splitext(sample)[0]

### init graph
graph = Graph(directed=True)


### set graph properties
for (name, prop, prop_type) in noderad_graph_operations.get_new_properties("init-graph"):
    if name.startswith("v_"):
        vars()[name] = graph.new_vertex_property(prop_type)
        graph.vertex_properties[prop] = vars()[name]
    if name.startswith("e_"):
        vars()[name] = graph.new_edge_property(prop_type)
        graph.edge_properties[prop] = vars()[name]


### create first empty node for graph
node = graph.add_vertex()
v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=sample)
v_name[node] = ""
v_seq[node] = ""
v_qual[node] = ""
nodes = [node]
edges = []

### add vertices
def add_verticles(_reads, idx_nodes):
    for record in SeqIO.parse(_reads, format):
        idx_nodes += 1
        node = graph.add_vertex()
        v_id[node] = "{idx}_{sample}".format(idx=idx_nodes, sample=sample)
        v_name[node] = (record.id).split(' ', 1)[0]
        v_seq[node] = record.seq
        v_qual[node] = [10**(-1 * i / 10) for i in record.letter_annotations["phred_quality"]]
        v_q_qual[node] = record.letter_annotations["phred_quality"]
        nodes.append(node)

if reads.endswith(".gz"):
    with gzip.open(reads, "rt") as _reads:
        add_verticles(_reads, idx_nodes)
else:
    with open(reads, "rU") as _reads:
        add_verticles(_reads, idx_nodes)

### add edges from all-vs-all alignment of reads (please see rule minimap2)
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
            query_node = find_vertex(graph, v_name, read.query_name)[0]
            ref_node = find_vertex(graph, v_name, read.reference_name)[0]
            if (not graph.vertex_index[ref_node] in graph.get_all_neighbors(query_node)):
                # add edge from alignment in sam-file
                edge = graph.add_edge(query_node, ref_node)
                e_dist[edge] = nm
                e_cs[edge] = cig

                # sam format properties to write results of optimal solution to sam file (rule optimized_solution)
                e_sam_cigar[edge] = read.cigarstring
                e_sam_qname[edge] = read.query_name
                e_sam_flag[edge] = read.flag
                e_sam_rname[edge] = read.reference_id
                e_sam_pos[edge] = read.reference_start
                e_sam_mapq[edge] = read.mapping_quality
                e_sam_cigar[edge] = read.cigar
                e_sam_rnext[edge] = read.next_reference_id
                e_sam_pnext[edge] = read.next_reference_start
                e_sam_tlen[edge] = read.template_length
                e_sam_seq[edge] = read.query_sequence
                e_sam_all_tags[edge] = read.tags

                qual_idx = 0
                for (op, length) in read.cigartuples:
                    mutrate = mut_total
                    qual_query_node = graph.vertex_properties["quality"][query_node]
                    if op == 7 or op == 0:  # on match
                        for i in qual_query_node[qual_idx:length]:
                            likelihood *= 1 - i
                        qual_idx = length
                    if op == 8 or op == 1 or op == 2 or op == 4:
                        if op == 8 or op == 4:  # on mismatch: substitution/snp or on softclip
                            if mut_subst:
                                mutrate = mut_subst
                                n_snp += length
                        if op == 1:  # on mismatch: insertion
                            if mut_ins:
                                mutrate = mut_ins
                            n_ins += length
                        if op == 2:  # on mismatch: deletion
                            if mut_del:
                                mutrate = mut_del
                            n_del += length
                        for i in qual_query_node[qual_idx:length]:
                            likelihood *= (1 - mutrate) * float(1/3) * i + mutrate * (1 - i)
                        qual_idx = length
                e_lh[edge] = likelihood
                edges.append((graph.vertex_index[edge.source()], graph.vertex_index[edge.target()]))
sam.close()

graph.remove_vertex(0)
n_nodes = len(nodes) - 1
n_edges = len(edges)
sys.stderr.write("graph construction summary for sample {}:\n nodes:\t{}\n edges:\t{}\n max distance:\t{}".format(sample, n_nodes, n_edges, max_NM))
sys.stderr.write("\n\tmutations:\n\t\tsnp's and substitutions: {}\n\t\tinsertions: {}\n\t\tdeletions: {}\n\n".format(n_snp, n_ins, n_del))

graph.save(snakemake.output.get("graph_xml"))
pos = sfdp_layout(graph)
if output_figure:
    graph_draw(graph, vertex_color=[1, 1, 1, 0],
               edge_color=graph.edge_properties["likelihood"],
               pos=pos, vertex_size=1, output=output_figure)
