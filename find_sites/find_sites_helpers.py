import math
import re
import time

from Bio import Phylo
import pandas as pd

import config
import utils


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

def get_seed_window_dict(seeds):
    seed_window_dict = {}
    for seed in seeds:
        site = utils.reverse_complement(seed)
        window_dict = {}
        window_dict[site+'A'] = '8mer-1a'
        for nt in ['C','U','G','X']:
            window_dict[site+nt] = '7mer-m8'
        other_nts = ['A','U','C','G','X']
        other_nts.remove(site[0])
        for nt1 in other_nts:
            window_dict[nt1+site[1:]+'A'] = '7mer-1a'
            for nt8 in ['C','U','G','X']:
                window_dict[nt1+site[1:]+nt8] = '6mer'
        seed_window_dict[seed] = window_dict
    return seed_window_dict

def get_matches(locs1,locs2,sp):
    l1,l2 = 0,0
    return_list = [None]*len(locs1)
    while (l1 < len(locs1)) & (l2 < len(locs2)):
        diff = (locs1[l1] - locs2[l2])
        if diff > 1:
            l2 += 1
        elif diff < -1:
            l1 += 1
        else:
            return_list[l1] = sp
            l1 += 1
            l2 += 1
    return return_list


def find_aligning_species(utr_df,seedm8,species,num_sites):
    if (len(species) == 0) | (config.REF_SPECIES not in species):
        return ['']*num_sites

    utr_df = utr_df.groupby('Species').first()
    ref_utr = utr_df.loc[config.REF_SPECIES]['UTR sequence']

    site = utils.reverse_complement(seedm8)
    regex = '-*'.join(site[1:])

    # ref_sites = [(m.start(),m.start() + re.match(regex,ref_utr[m.start():]).end()-1,[]) for m in re.finditer('(?={})'.format(regex), ref_utr)]
    ref_sites = [m.start() for m in re.finditer('(?={})'.format(regex), ref_utr)]

    species_list = []

    # loop through species that have this miRNA
    for sp in species:

        # species is not in alignment
        if sp not in utr_df.index:
            continue

        # species is the reference species, don't need to recalculate sites -> we know they all align to the reference
        elif sp == config.REF_SPECIES:
            species_list.append([sp]*len(ref_sites))

        # look for sites in this species that align to sites in the reference species
        else:
            sites = [m.start() for m in re.finditer('(?={})'.format(regex), utr_df.loc[sp]['UTR sequence'])]
            if len(sites) != 0:
                species_list.append(get_matches(ref_sites,sites,sp))


    if len(species_list) == 0:
        return ['']*num_sites


    species_list = zip(*species_list)
    species_list = [[x for x in sp if x is not None] for sp in species_list]
    # species_list = [';'.join(sp) for sp in species_list]
    return species_list
    

def get_site_info(utr_no_gaps,seed,window_dict):
    # print utr_no_gaps
    # utr_no_gaps = 'X'+utr_no_gaps.replace('\n','')+'X'
    rc_seed = utils.reverse_complement(seed)
    locs = [m.start() for m in re.finditer('(?={})'.format(rc_seed[1:]),utr_no_gaps)]

    windows = [utr_no_gaps[x-1:x+7] for x in locs]
    site_types = [window_dict[x] for x in windows]

    site_starts = [x - 1 - int(y in ['8mer-1a','7mer-m8']) for (x,y) in zip(locs,site_types)]
    site_ends = [x +5+ int(y in ['8mer-1a','7mer-1a']) for (x,y) in zip(locs,site_types)]
    
    return tuple(site_starts),tuple(site_ends),tuple(site_types)


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

# def get_branch_length_score(species_list,bin):
#     species_list = [x for x in species_list if x != -1]
#     if len(species_list) == 0:
#         return 0
#     all_paths = []
#     for sp in species_list:
#         all_paths += config.TREES[bin][sp]
#     all_paths = list(set(list(all_paths)))
#     bls = sum([x[1] for x in all_paths])
#     return bls

# def calculate_pct(aligning_species,bin,site_type,seed):
#     if site_type == '6mer':
#         return 0.0,0

#     bls = get_branch_length_score(aligning_species,bin)
#     b0,b1,b2,b3 = config.PARAMS[site_type][seed]

#     score = max(0.0, b0 + b1/(1.0 + (math.e*(0.0-b2)*bls) + b3))

#     conserved = int(score >= config.CONSERVATION_CUTOFFS[site_type])

#     return [score, conserved]



# seed = 'GAGGUAG'
# seed_window_dict = get_seed_window_dict([seed])
# #      0         1         2         3         4         5         6         7   
# #      xxxxx-xx                   o-ooooo           ooooooo                     ooooooo     
# utr = 'CUACC-UCAAUUUGGA-CUCCCUUUCAU-ACCUCUGGAA-GAGCUUACCUCAUCCUAGCAUAAAAUCUUUGCACUACCUC'
# utr2 = 'CUACG-UCAAUUUGGA-CUCCCUUUCAU-UCCUCUGGAA-GAGCUUACCUCAUCCUAGCAUAAAAUCUUUGCACUAGCUC'
# utr_df = pd.DataFrame({'Species': ['9606','9615'],'UTR sequence': [utr,utr2]})
# utr_df = utr_df.set_index('Species')
# utr_no_gaps = utr.replace('-','')
# print find_sites(utr_df,seed)
# window_dict = seed_window_dict[seed]


