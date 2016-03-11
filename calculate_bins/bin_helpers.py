import bisect
import math
import sys
import time

from Bio import Phylo
import numpy as np

import config
from utils import reverse_complement


def parse_tree(tree_file, ref_species):
    """
    Parses a phylogenetic tree and find the edge path between the reference
    species and all the other species in the tree

    Parameters
    ----------
    tree_file: string, name of file that contains the phylogenetic tree
        tree must be in newick format

    ref_species: string, species id of the reference species
        must be contained in the tree

    Output
    ------
    dictionary: for every species, gives the edge path
        to the reference species
    """
    # read in the phylogenetic tree using Phylo
    tree = Phylo.read(tree_file, 'newick')

    # gene all the leaves of the tree and check for the reference species
    leaves = tree.get_terminals()
    for i, leaf in enumerate(leaves):
        if leaf.name == ref_species:
            ref_leaf = leaves.pop(i)
    if ref_leaf.name != ref_species:
        print 'Error: reference species not in tree'
        sys.exit()

    # dictionary that stores the paths to the reference species
    species_to_path = {ref_leaf.name: [(ref_leaf.name, 0)]}

    # give each node a number so that we can uniquely identify edges
    for i, node in enumerate(tree.get_nonterminals()):
        node.name = 'n{}'.format(i)

    # for all the species, find the path to the reference species
    for leaf in leaves:
        path = tree.trace(ref_leaf, leaf)
        lca = tree.common_ancestor(ref_leaf, leaf).name
        path_lengths = [(ref_leaf.name, float(ref_leaf.branch_length))]
        for edge in path:
            if edge.name != lca:
                path_lengths.append((edge.name, float(edge.branch_length)))

        # store path information in the dictionary
        species_to_path[leaf.name] = path_lengths

    return species_to_path


def compare_seqs(seq, ref_utr, species_id):
    """
    Compares two aligned sequences and find where they match

    Parameters
    ----------
    seq: string, sequence to compare

    ref_utr: string, sequence of reference species

    species id: string, species id of the sequence to be compared

    Output
    ------
    list: for every position in the reference utr, list the species id
        if the sequences match, otherwise list -1
    """

    # convert the sequence to upper case
    # this should already be done for the reference utr
    seq = seq.upper()

    # truncate or add gaps to the sequence so that it matches the
    # reference utr in length
    ref_length = len(ref_utr)
    if len(seq) > ref_length:
        seq = seq[:ref_length]
    elif len(seq) < ref_length:
        seq = seq + '-'*(ref_length-len(seq))

    return [species_id if int(x == y) else -1 for (x, y) in zip(seq, ref_utr)]


def get_branch_length_score(species_list):
    """
    Given a list of species, compute the branch length score

    Parameters
    ----------
    species_list: list of strings, contains species ids

    Output
    ------
    float: branch length score for this list of species
    """

    # species list may contain some -1, remove these first
    species_list = [x for x in species_list if x != -1]

    # if the species list is empty, return a bls of 0
    if len(species_list) == 0:
        return 0

    # get all the edges that connect these species to the reference
    all_paths = []
    for sp in species_list:
        all_paths += config.SPECIES_TO_PATH[sp]

    # get a unique list of these edges and compute bls
    all_paths = list(set(list(all_paths)))
    bls = sum([x[1] for x in all_paths])

    return bls


def get_median_binned_list(counts, values, total, num_gap):
    """
    Given a list of values and counts, return the median without gaps
    Example: original list = [1, 1, -, 2, 3, -, 3, 3, 4]
        input into function:
            values = [0, 1, 2, 3, 4]
            counts = [2, 2, 1, 3, 1]
            total = 7
            num_gap = 2 (listed as 0 in the values list)
        output:
            median = 3

    Parameters
    ----------
    counts: list of ints

    values: list of floats, sorted

    total: total number of values

    num_gap: number of gaps in original list

    Output
    ------
    float: median of original list
    """
    # if the total number of values is even, we have to average 2 values
    if total % 2 == 0:
        running = 0

        # iterate through values until we have accumulated half the list
        for i, c in enumerate(counts):
            running += c
            if running > (total / 2) + num_gap:
                return values[i]
            if running == (total / 2) + num_gap:
                return (values[i] + values[i + 1]) / 2.0

    # otherwise, just find the middle value
    else:
        running = 0
        for i, c in enumerate(counts):
            running += c
            if running >= math.ceil(total / 2.0) + num_gap:
                return values[i]

    # return 0 if we cannot find the median for some reason
    return 0


def get_bin(bls):
    """
    Given a bls, assign to a bin using precomputed thresholds

    Parameters
    ----------
    bls: float

    Output
    ------
    int: bin
    """

    # get thresholds from config file
    thresholds = config.BL_THRESHOLDS

    # use binary search to find biggest threshold <= to the bls
    return bisect.bisect_right(thresholds, bls)
