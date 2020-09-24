import sys
import mappy as mp  # Minimap2 library (Python binding), see https://github.com/lh3/minimap2/tree/master/python

sys.stderr = open(snakemake.log[0], "w")
alignment = mp.Aligner(snakemake.input[0])  # load or build index
if not alignment:
    raise Exception("ERROR: Failed to build alignment index.")

n_nodes = len(alignment.seq_names)  # number of nodes
n_edges = 0  # number of edges # number of edges
max_NM = 0  # highest edit distance
thread_buffer = mp.ThreadBuffer()

for name, seq, qual in mp.fastx_read(snakemake.input[0]): # read a fasta/q sequence
    for hit in alignment.map(seq, buf=thread_buffer): # traverse alignments
        if not name == hit.ctg:
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
x = list(mp.fastx_read(snakemake.input[0]))[3]  ### indexing on type: generator ###
print(x) ###

with open(snakemake.output[0], "w") as out:  ###
    print("Success!", file=out)  ###

########################################################################################
