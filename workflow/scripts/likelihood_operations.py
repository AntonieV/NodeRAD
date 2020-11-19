from graph_tool.all import *
from itertools import zip_longest, combinations_with_replacement
from collections import Counter
import math
import editdistance


# if not reverse: calculates likelihood for an edge (from ref to query) from the cigartuples and quality of the query node
# if reverse: calculates likelihood for backwards an edge (from ref to query) from the cigartuples and ref quality of an existing edge (from query to ref)
def get_alignment_likelihood(cigar_tuples, quality, mut_rates, stats=None, reverse=False):
    qual_idx = 0
    likelihood = 1.0
    for (op, length) in cigar_tuples:
        mutrate = mut_rates["mut_total"]
        if op == 7 or op == 0:  # on match
            for i in quality[qual_idx:length]:
                likelihood *= 1 - i
            qual_idx = length
        if op == 8 or op == 1 or op == 2 or op == 4:
            if op == 8 or op == 4:  # on mismatch: substitution/snp or on softclip
                if mut_rates["mut_subst"]:
                    mutrate = mut_rates["mut_subst"]
                if stats:
                    stats["n_snp"] += length
            if op == 1:
                if mut_rates["mut_ins"]:
                    # on mismatch of not reversed edge: insertion in query node
                    # on mismatch of backwards edge: insertion in reference node corresponds to deletion in query node
                    mutrate = mut_rates["mut_ins"] if not reverse else mut_rates["mut_del"]
                if stats:
                    stats["n_ins"] += length
            if op == 2:
                if mut_rates["mut_del"]:
                    # on mismatch of not reversed edge: deletion in query node
                    # on mismatch of backwards edge: deletion in reference node corresponds to insertion in query node
                    mutrate = mut_rates["mut_del"] if not reverse else mut_rates["mut_ins"]
                if stats:
                    stats["n_del"] += length
            for i in quality[qual_idx:length]:
                likelihood *= (1 - mutrate) * float(1 / 3) * i + mutrate * (1 - i)
            qual_idx = length
    if stats:
        return likelihood, stats
    return likelihood


def get_cigar_tuples(comp, node, allele):
    source_nodes = [comp.vertex_index[node] for node in
                    find_vertex(comp, comp.vp["sequence"], comp.vp["sequence"][node])]
    target_nodes = [comp.vertex_index[node] for node in find_vertex(comp, comp.vp["sequence"], allele)]
    for node_s in source_nodes:
        for node_t in target_nodes:
            edge = comp.edge(node_s, node_t)
            if edge:
                return comp.edge_properties["cigar-tuples"][edge], False
            rev_edge = comp.edge(node_t, node_s)
            if rev_edge:
                return comp.edge_properties["cigar-tuples"][rev_edge], True
    return 0


# returns a map with lists of all possible combinations of fractions from a given number of alleles n and
# the ploidy given in config file
def get_candidate_vafs(n, ploidy):
    n_alleles = n + n % ploidy  # each locus is of same ploidy
    get_combination_vafs = lambda comb: [sum(x == i for x in comb) / n_alleles for i in range(n)]
    return map(get_combination_vafs, combinations_with_replacement(range(n), n_alleles))


def get_candidate_alleles(graph, reads, noise):
    read_support = Counter()
    for read in reads:
        read_support[graph.vertex_properties['sequence'][read]] += 1
    return sorted(seq for seq, count in read_support.items() if count >= noise)


def get_allele_likelihood_read(allele, comp, node, mut_rates):
    # obtian one arbitrary out edge of node that points to another node with sequence = allele
    out_neighbors = comp.get_out_neighbors(node)
    in_neighbors = comp.get_in_neighbors(node)
    # search for existing edge from given node to target node with sequence = allele
    for neighbor in out_neighbors:
        if comp.vertex_properties['sequence'][neighbor] == allele:
            return comp.edge_properties["likelihood"][comp.edge(node, neighbor)]
    # search for reverse target node with sequence = allele to given node
    qual = comp.vp["quality"][node]
    for neighbor in in_neighbors:
        if comp.vertex_properties['sequence'][neighbor] == allele:
            return get_alignment_likelihood(list(eval(comp.edge_properties['cigar-tuples'][comp.edge(neighbor, node)])),
                                                qual, mut_rates, reverse=True)
    cigar_tuples = get_cigar_tuples(comp, node, allele)
    if cigar_tuples:
        return get_alignment_likelihood(list(eval(cigar_tuples[0])), qual, mut_rates, reverse=cigar_tuples[1])
    return 0


def calc_vafs_likelihood_read(vafs, comp, node, alleles, mut_rates):
    vafs_lh_read = sum(
        vaf * get_allele_likelihood_read(allele, comp, node, mut_rates)
        for vaf, allele in zip(vafs, alleles)
    )
    return math.log(vafs_lh_read) if vafs_lh_read > 0 else -float("inf")


def calc_vafs_likelihood(vafs, comp, nodes, alleles, mut_rates):
    return sum(calc_vafs_likelihood_read(vafs, comp, node, alleles, mut_rates) for node in nodes)


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def get_candidate_loci(n, ploidy):
    n_alleles = n + n % ploidy  # each locus is of same ploidy
    get_combination_loci = lambda comb: list(grouper(comb, ploidy))
    return map(get_combination_loci, combinations_with_replacement(range(n), n_alleles))


def get_alleles_edit_distance(comp, allele_1, allele_2):
    nodes_allele_1 = [comp.vertex_index[node] for node in
                    find_vertex(comp, comp.vp["sequence"], allele_1)]
    nodes_allele_2 = [comp.vertex_index[node] for node in find_vertex(comp, comp.vp["sequence"], allele_2)]
    for node_s in nodes_allele_1:
        for node_t in nodes_allele_2:
            if comp.edge(node_s, node_t):
                return comp.edge_properties["distance"][comp.edge(node_s, node_t)]
            if comp.edge(node_t, node_s):
                return comp.edge_properties["distance"][comp.edge(node_t, node_s)]
    return editdistance.eval(allele_1, allele_2)


def get_allele_likelihood_allele(comp, locus_alleles):
    # get subgraph induced by locus alleles (using only ONE node per allele)
    # get minimum spanning tree (with edit distance as weight)
    for allele_1 in locus_alleles:
        for allele_2 in locus_alleles:
            if not allele_1 == allele_2:
                dist = get_alleles_edit_distance(comp, allele_1, allele_2)
                print(dist)

    # spanning_tree = ...
    # likelihood = 0
    # for root in spanning_tree.nodes():
    #     rooted_mst_likelihood = 0
    #     # traverse spanning tree from root
    #     # for each traversed edge, calculate pairHMM probability from parent to child node
    #     # and sum up in log space
    #     rooted_mst_likelihood += ...
    # # add to overall sum
    # return likelihood


# def get_indicator_constrait():



def calc_loci_likelihoods(max_likelihood_vafs, comp, alleles, loci):
    # formula from slide 9
    indicator_constraint = get_indicator_constrait()
    if indicator_constraint:
        get_allele_likelihood_allele(comp, alleles)
        # calculate double product from slide 9
        return ...
    else:
        return 0

