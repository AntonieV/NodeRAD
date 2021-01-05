from graph_tool.all import *
from itertools import zip_longest, combinations_with_replacement, chain, combinations, permutations
from collections import Counter
import math


# if not reverse: calculates likelihood for an edge (from ref to query) from the cigartuples and quality of the query node
# if reverse: calculates likelihood for backwards an edge (from ref to query) from the cigartuples and ref quality of an existing edge (from query to ref)
def get_alignment_likelihood(graph, cigar_tuples, qual, reverse=False):
    quality = [10 ** (-1 * i / 10) for i in qual]
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
            for i in quality[qual_idx:(qual_idx+length)]:
                likelihood *= 1 - i
            qual_idx += length
        if op == 8 or op == 1 or op == 2 or op == 4:
            mut_rate = 0
            if op == 8 or op == 4:  # on substitution/snp or on softclip mismatch in cigar-string
                mut_rate = _subst
            if op == 1:  # on insertion mismatch in cigar-string
                mut_rate = _ins
            if op == 2:  # on deletion mismatch in cigar-string
                mut_rate = _del
            for i in quality[qual_idx:(qual_idx+length)]:
                likelihood *= (1 - mut_rate) * float(1 / 3) * i + mut_rate * (1 - i)
            qual_idx += length
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
                heterozyg_factor = _subst
            if op == 1:  # on insertion mismatch in cigar-string
                heterozyg_factor = _ins
            if op == 2:  # on deletion mismatch in cigar-string
                heterozyg_factor = _del
            heterozygosity *= heterozyg_factor ** length
    return heterozygosity


def get_cigar_tuples(comp, seq, allele):
    source_nodes = [comp.vertex_index[node] for node in
                    find_vertex(comp, comp.vp["sequence"], seq)]
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


# returns a map with lists of all possible combinations of fractions from a given number of alleles n and
# the ploidy given in config file
def get_candidate_vafs(n, ploidy):
    n_alleles = get_max_parsimony_n_alleles(n, ploidy)
    get_combination_vafs = lambda comb: [sum(x == i for x in comb) / n_alleles for i in range(n)]
    return map(get_combination_vafs, combinations_with_replacement(range(n), n_alleles))


def get_candidate_alleles(comp, reads, noise, cluster_size):
    read_support = Counter()
    for read in reads:
        read_support[comp.vertex_properties['sequence'][read]] += 1
    # filters noise for clusters above a certain cluster size
    # treshold value of noise and cluster-size are given in config
    if comp.num_vertices() >= cluster_size:
        return sorted(seq for seq, count in read_support.items() if count >= noise)
    return sorted(seq for seq, count in read_support.items())


def get_allele_likelihood_read(comp, allele, node, read_allele_likelihoods):
    # obtian one arbitrary out edge of node that points to another node with sequence = allele
    out_neighbors = comp.get_out_neighbors(node)
    qual = comp.vp["quality-q-vals"][node]
    # search for existing edge from given node to target node with sequence = allele
    for neighbor in out_neighbors:
        seq_neighbor = comp.vertex_properties['sequence'][neighbor]
        if seq_neighbor == allele:
            if not (seq_neighbor, allele) in read_allele_likelihoods:
                read_allele_likelihoods[(seq_neighbor, allele)] = get_alignment_likelihood(comp, list(eval(comp.edge_properties['cigar-tuples'][comp.edge(node, neighbor)])), qual, reverse=False)
            return read_allele_likelihoods[(seq_neighbor, allele)]
    # search for reverse target node with sequence = allele to given node
    in_neighbors = comp.get_in_neighbors(node)
    for neighbor in in_neighbors:
        seq_neighbor = comp.vertex_properties['sequence'][neighbor]
        if seq_neighbor == allele:
            if not (seq_neighbor, allele) in read_allele_likelihoods:
                read_allele_likelihoods[(seq_neighbor, allele)] = get_alignment_likelihood(comp, list(eval(comp.edge_properties['cigar-tuples'][comp.edge(neighbor, node)])), qual, reverse=True)
            return read_allele_likelihoods[(seq_neighbor, allele)]
    return 0


def calc_vafs_likelihood_read(comp, vafs, node, alleles, read_allele_likelihoods):
    vafs_lh_read = sum(
        vaf * get_allele_likelihood_read(comp, allele, node, read_allele_likelihoods)
        for vaf, allele in zip(vafs, alleles)
    )
    return math.log(vafs_lh_read) if vafs_lh_read > 0 else -float("inf")


def calc_vafs_likelihood(comp, vafs, nodes, alleles, read_allele_likelihoods):
    return sum(calc_vafs_likelihood_read(comp, vafs, node, alleles, read_allele_likelihoods) for node in nodes)


# source: https://docs.python.org/3/library/itertools.html
# Collect data into fixed-length chunks or blocks"
# grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def get_candidate_loci(n, ploidy):
    n_alleles = get_max_parsimony_n_alleles(n, ploidy)
    get_combination_loci = lambda comb: tuple(sorted(map(lambda l: tuple(sorted(l)), grouper(comb, ploidy))))
    return set(get_combination_loci(p) for comb in combinations_with_replacement(range(n), n_alleles) for p in
               permutations(comb))


def get_allele_likelihood_allele(comp, loci_alleles):
    likelihood = 1.0
    for locus_alleles in loci_alleles:
        # heterozygosity for all combinations of allele pairs
        for allele_i, allele_j in list(combinations(locus_alleles, 2)):
            cigar = get_cigar_tuples(comp, allele_i, allele_j)
            if cigar:
                likelihood += math.log(get_heterozygosity(comp, list(eval(cigar[0])), reverse=cigar[1]))
    return math.exp(likelihood)


def indicator_constrait(ploidy, max_likelihood_vafs, loci):
    count_alleles = Counter(list(chain(*loci)))
    return all([count_alleles[idx] == max_likelihood_vafs[idx] * len(loci) * ploidy for locus in loci for idx in locus])


def calc_loci_likelihoods(comp, max_likelihood_vafs, alleles, loci):
    ploidy = comp.gp["ploidy"]
    if indicator_constrait(ploidy, max_likelihood_vafs, loci):
        loci_alleles = list(list(map(lambda x: alleles[x], locus)) for locus in loci)
        return get_allele_likelihood_allele(comp, loci_alleles)
    return 0


def get_sorted_loci_alleles(alleles, max_likelihood_loci):
    return sorted(set([alleles[loc] for loci in max_likelihood_loci for loc in loci]))


def get_alleles_matched_to_loci(alleles, max_likelihood_loci):
    matched_al_loc = []
    for locus in max_likelihood_loci:
        matched_al_loc.append(tuple(alleles[idx] for idx in locus))
    return matched_al_loc


def get_gt_indices(alleles, max_likelihood_loci, loc_alleles):
    matched_al_loc = get_alleles_matched_to_loci(alleles, max_likelihood_loci)
    gt_indices = []
    for locus in matched_al_loc:
        gt_indices.append(tuple(loc_alleles.index(allele) for allele in locus))
    return gt_indices


def get_genotype(locus):
    return '/'.join(map(str, [idx for idx in locus]))

