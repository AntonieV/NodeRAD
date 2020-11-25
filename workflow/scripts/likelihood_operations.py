from graph_tool.all import *
from itertools import zip_longest, combinations_with_replacement
from collections import Counter
import math
import graph_operations


# if not reverse: calculates likelihood for an edge (from ref to query) from the cigartuples and quality of the query node
# if reverse: calculates likelihood for backwards an edge (from ref to query) from the cigartuples and ref quality of an existing edge (from query to ref)
def get_alignment_likelihood(graph, cigar_tuples, quality, stats=None, reverse=False):
    qual_idx = 0
    likelihood = 1.0
    _subst = graph.gp["subst-mutation-rates"]
    # mismatch of not reversed edge: insertion in query node
    _ins = graph.gp["ins-mutation-rates"]
    # mismatch of not reversed edge: deletion in query node
    _del = graph.gp["del-mutation-rates"]
    if reverse:
        # on mismatch of backwards edge: deletion in reference node corresponds to insertion in query node
        _ins = _del
        # on mismatch of backwards edge: insertion in reference node corresponds to deletion in query node
        _del = graph.gp["ins-mutation-rates"]

    for (op, length) in cigar_tuples:
        if op == 7 or op == 0:  # on match in cigar-string
            for i in quality[qual_idx:length]:
                likelihood *= 1 - i
            qual_idx = length
        if op == 8 or op == 1 or op == 2 or op == 4:
            mut_rate = 0
            if op == 8 or op == 4:  # on substitution/snp or on softclip mismatch in cigar-string
                mut_rate = _subst
                if stats:
                    stats["n_snp"] += length
            if op == 1:  # on insertion mismatch in cigar-string
                mut_rate = _ins
                if stats:
                    stats["n_ins"] += length
            if op == 2:  # on deletion mismatch in cigar-string
                mut_rate = _del
                if stats:
                    stats["n_del"] += length
            for i in quality[qual_idx:length]:
                likelihood *= (1 - mut_rate) * float(1 / 3) * i + mut_rate * (1 - i)
            qual_idx = length
    if stats:
        return likelihood, stats
    return likelihood


def get_heterozygosity(comp, cigar_tuples, reverse=False):
    heterozygosity = 1.0
    heterozyg_factor = 1.0
    _subst = comp.gp["subst-heterozygosity"]
    # mismatch of not reversed edge: insertion in query node
    _ins = comp.gp["ins-heterozygosity"]
    # mismatch of not reversed edge: deletion in query node
    _del = comp.gp["del-heterozygosity"]
    if reverse:
        # on mismatch of backwards edge: deletion in reference node corresponds to insertion in query node
        _ins = _del
        # on mismatch of backwards edge: insertion in reference node corresponds to deletion in query node
        _del = comp.gp["ins-heterozygosity"]

    heterozyg_all = _subst + _ins + _del

    for (op, length) in cigar_tuples:
        if op == 7 or op == 0:  # on match in cigar-string
            heterozygosity *= (1 - heterozyg_all) ** length
        if op == 8 or op == 1 or op == 2 or op == 4:
            if op == 8 or op == 4:  # on substitution/snp or on softclip mismatch in cigar-string
                heterozyg_factor = _subst ** length
            if op == 1:  # on insertion mismatch in cigar-string
                heterozyg_factor = _ins
            if op == 2:  # on deletion mismatch in cigar-string
                heterozyg_factor = _del
            heterozygosity *= heterozyg_factor
    return heterozygosity


def get_cigar_tuples(comp, seq_node, allele):
    source_nodes = [comp.vertex_index[node] for node in
                    find_vertex(comp, comp.vp["sequence"], seq_node)]
    target_nodes = [comp.vertex_index[node] for node in find_vertex(comp, comp.vp["sequence"], allele)]
    for node_s in source_nodes:
        for node_t in target_nodes:
            edge = comp.edge(node_s, node_t)
            if edge:
                return comp.edge_properties["cigar-tuples"][edge], False
            rev_edge = comp.edge(node_t, node_s)
            if rev_edge:
                return comp.edge_properties["cigar-tuples"][rev_edge], True
    return None


def get_max_parsimony_n_alleles(n, ploidy):
    """Get the expected number of true alleles in the connected component
    given that n alleles have been observed.

    The returned value can be seen as the parsimonious solution: if ploidy is
    higher than the number of observed alleles, there must be at least as many
    alleles as the ploidy (though some are the same, e.g. at homozygous loci).
    If ploidy is less than the observed alleles, there must be so many additional
    alleles that are the same as one of the observed such that the number of alleles
    can be divided by the ploidy."""
    remainder = n % ploidy
    if ploidy >= n:
        return ploidy
    elif remainder:
        return n + ploidy - remainder
    else:
        return n

    # if ploidy >= n:
    #     return ploidy
    # elif (remainder := n % ploidy):
    #     return n + ploidy - remainder
    # else:
    #     return n


# returns a map with lists of all possible combinations of fractions from a given number of alleles n and
# the ploidy given in config file
def get_candidate_vafs(n, ploidy):
    n_alleles = get_max_parsimony_n_alleles(n, ploidy)
    get_combination_vafs = lambda comb: [sum(x == i for x in comb) / n_alleles for i in range(n)]
    return map(get_combination_vafs, combinations_with_replacement(range(n), n_alleles))


def get_candidate_alleles(graph, reads, noise):
    read_support = Counter()
    for read in reads:
        read_support[graph.vertex_properties['sequence'][read]] += 1
    return sorted(seq for seq, count in read_support.items() if count >= noise)


def get_allele_likelihood_read(comp, allele, node):
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
            return get_alignment_likelihood(comp, list(eval(comp.edge_properties['cigar-tuples'][comp.edge(neighbor, node)])),
                                                qual, reverse=True)
    cigar_tuples = get_cigar_tuples(comp, comp.vp["sequence"][node], allele)
    if cigar_tuples:
        return get_alignment_likelihood(comp, list(eval(cigar_tuples[0])), qual, reverse=cigar_tuples[1])
    return 0


def calc_vafs_likelihood_read(comp, vafs, node, alleles):
    vafs_lh_read = sum(
        vaf * get_allele_likelihood_read(comp, allele, node)
        for vaf, allele in zip(vafs, alleles)
    )
    return math.log(vafs_lh_read) if vafs_lh_read > 0 else -float("inf")


def calc_vafs_likelihood(comp, vafs, nodes, alleles):
    return sum(calc_vafs_likelihood_read(comp, vafs, node, alleles) for node in nodes)


# source: https://docs.python.org/3/library/itertools.html
def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def get_candidate_loci(n, ploidy):
    n_alleles = get_max_parsimony_n_alleles(n, ploidy)
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
    return None


# get subgraph induced by locus alleles (using only ONE node per allele)
def init_alleles_subgraph():
    allele_subgraph = Graph(directed=False)
    # set graph properties
    for (name, prop, prop_type) in graph_operations.set_properties("alleles"):
        if name.startswith("g_"):
            vars()[name] = allele_subgraph.new_graph_property(prop_type)
            allele_subgraph.graph_properties[prop] = vars()[name]
        if name.startswith("v_"):
            vars()[name] = allele_subgraph.new_vertex_property(prop_type)
            allele_subgraph.vertex_properties[prop] = vars()[name]
        if name.startswith("e_"):
            vars()[name] = allele_subgraph.new_edge_property(prop_type)
            allele_subgraph.edge_properties[prop] = vars()[name]
    return allele_subgraph


# get minimum spanning tree (with edit distance as weight)
def get_allele_subgraph(comp, alleles, dir_subgraph, comp_nr):
    allele_subgraph = init_alleles_subgraph()
    allele_subgraph.graph_properties["subst-heterozygosity"] = comp.graph_properties["subst-heterozygosity"]
    allele_subgraph.graph_properties["ins-heterozygosity"] = comp.graph_properties["ins-heterozygosity"]
    allele_subgraph.graph_properties["del-heterozygosity"] = comp.graph_properties["del-heterozygosity"]
    for allele_1 in alleles:
        node_exists = find_vertex(allele_subgraph, allele_subgraph.vp["allele-sequence"], allele_1)
        if node_exists:
            node_1 = node_exists[0]
        else:
            node_1 = allele_subgraph.add_vertex()
            allele_subgraph.vp["allele-sequence"][node_1] = allele_1

        for allele_2 in alleles:
            if not allele_1 == allele_2:
                node_exists = find_vertex(allele_subgraph, allele_subgraph.vp["allele-sequence"], allele_2)
                if node_exists:
                    node_2 = node_exists[0]
                else:
                    node_2 = allele_subgraph.add_vertex()
                    allele_subgraph.vp["allele-sequence"][node_2] = allele_2
                edit_distances = get_alleles_edit_distance(comp, allele_1, allele_2)
                if edit_distances:
                    edge = allele_subgraph.add_edge(node_1, node_2)
                    allele_subgraph.ep["allele-distance"][edge] = edit_distances

    if dir_subgraph:
        graph_operations.set_dir(dir_subgraph)
        xml = "{}/{}.subgraph.xml.gz".format(dir_subgraph, str(comp_nr))
        figure = "{}/{}.subgraph.pdf".format(dir_subgraph, str(comp_nr))
        graph_operations.save_and_draw_graph(allele_subgraph,
                                             xml_out=xml,
                                             figure_out=figure,
                                             e_color="allele-distance",
                                             v_size=9)
    return allele_subgraph


def get_spanning_tree(allele_subgraph, dir_mst, comp_nr):
    spanning_tree_map = min_spanning_tree(allele_subgraph, weights=allele_subgraph.ep["allele-distance"])
    spanning_tree_view = GraphView(allele_subgraph, efilt=spanning_tree_map)
    if dir_mst:
        graph_operations.set_dir(dir_mst)
        xml = "{}/{}.spanning-tree.xml.gz".format(dir_mst, str(comp_nr))
        figure = "{}/{}.spanning-tree.pdf".format(dir_mst, str(comp_nr))
        graph_operations.save_and_draw_graph(spanning_tree_view,
                                             xml_out=xml,
                                             figure_out=figure,
                                             e_color="allele-distance",
                                             v_size=9)
    return Graph(spanning_tree_view, prune=True, directed=False)

def get_allele_likelihood_allele(comp, allele_subgraph, spanning_tree, locus_alleles):
    # traverse spanning tree from root
    # for each traversed edge, calculate pairHMM probability from parent to child node
    # and sum up in log space
    likelihood = 0
    for allele_i in locus_alleles:
        node_source = find_vertex(spanning_tree, spanning_tree.vp["allele-sequence"], allele_i)[0]
        for allele_j in locus_alleles:
            if not allele_i == allele_j:
                node_target = find_vertex(spanning_tree, spanning_tree.vp["allele-sequence"], allele_j)[0]
                nodes, edges = shortest_path(spanning_tree, source=node_source, target=node_target)
                rooted_mst_likelihood = 0
                for edge in edges:
                    allele_parent = spanning_tree.vertex_properties["allele-sequence"][edge.source()]
                    allele_child = spanning_tree.vertex_properties["allele-sequence"][edge.target()]
                    cigar = get_cigar_tuples(comp, allele_parent, allele_child)
                    rooted_mst_likelihood += math.log(get_heterozygosity(allele_subgraph, list(eval(cigar[0])), reverse=cigar[1]))
                # add to overall sum
                likelihood += rooted_mst_likelihood
    return likelihood


def indicator_constrait(max_likelihood_vafs, loci):
    return sum(max_likelihood_vafs[idx] for locus in loci for idx in locus) == 1


def calc_loci_likelihoods(comp, max_likelihood_vafs, alleles, loci, allele_subgraph, spanning_tree):
    if indicator_constrait(max_likelihood_vafs, loci):
        locus_alleles = [alleles[idx] for locus in loci for idx in locus]
        return get_allele_likelihood_allele(comp, allele_subgraph, spanning_tree, locus_alleles)
    return -float("inf")


