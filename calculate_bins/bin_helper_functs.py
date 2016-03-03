import bisect
import math
import time

from Bio import Phylo
import numpy as np

import config
from utils import reverse_complement


def parse_tree(tree_file,ref_species):
    tree = Phylo.read(tree_file,'newick')

    leaves = tree.get_terminals()
    for i,leaf in enumerate(leaves):
        if leaf.name == ref_species:
            ref_leaf = leaves.pop(i)
    assert ref_leaf.name == ref_species, "Error: reference species not in tree"
    species_to_path = {ref_leaf.name:[(ref_leaf.name,0)]}

    for i,node in enumerate(tree.get_nonterminals()):
        node.name = 'n{}'.format(i)

    for leaf in leaves:
        path = tree.trace(ref_leaf,leaf)
        lca = tree.common_ancestor(ref_leaf,leaf).name
        path_lengths = [(ref_leaf.name,float(ref_leaf.branch_length))]
        for edge in path:
            if edge.name != lca:
                path_lengths.append((edge.name,float(edge.branch_length)))
        species_to_path[leaf.name] = path_lengths
    return species_to_path

def get_median_binned_list(counts,values,total,num_gap):
    if total % 2 == 0:
        running = 0
        for i,c in enumerate(counts):
            running += c
            if running > (total/2)+num_gap:
                return values[i]
            if running == (total/2)+num_gap:
                return (values[i]+values[i+1])/2.0
        assert False, 'blargh'
        return 0
    else:
        running = 0
        for i,c in enumerate(counts):
            running += c
            if running >= math.ceil(total/2.0)+num_gap:
                return values[i]
        assert False, 'blargh'
        return 0

def compare_seqs(seq,ref_utr,species_id):
    seq = seq.upper()
    ref_length = len(ref_utr)
    if len(seq) > ref_length:
        seq = seq[:ref_length]
    elif len(seq) < ref_length:
        seq = seq + '-'*(ref_length-len(seq))
    return [species_id if int(x==y) else -1 for (x,y) in zip(seq,ref_utr)]

def get_branch_length_score(species_list):
    species_list = [x for x in species_list if x != -1]
    if len(species_list) == 0:
        return 0
    all_paths = []
    for sp in species_list:
        all_paths += config.SPECIES_TO_PATH[sp]
    all_paths = list(set(list(all_paths)))
    bls = sum([x[1] for x in all_paths])
    return bls

def get_bin(bls):
    '''Find smallest item greater-than or equal to bls.
    Raise ValueError if no such item exists.
    If multiple keys are equal, return the leftmost.

    '''
    thresholds = config.BL_THRESHOLDS
    i = bisect.bisect_left(thresholds, bls)
    if i == len(thresholds):
        return len(thresholds)
    if thresholds[i] == bls:
        return i+1
    return i
