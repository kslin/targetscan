import re
from string import maketrans
import sys

from Bio import Phylo
import numpy as np
import pandas as pd

import config


def rev_comp(seq):
    """
    Parameters:
    ==========
    seq: string, sequence to get reverse complement of

    Returns:
    =======
    float: reverse complement of seq in all caps
    """
    seq = seq.upper()
    intab = "AUCG"
    outtab = "UAGC"
    trantab = maketrans(intab, outtab)
    seq = seq[::-1]
    seq = seq.translate(trantab)
    return seq


def occurrences(string, sub):
    """
    Parameters:
    ==========
    string: string, longer string
    sub: string, substring you are looking for

    Returns:
    =======
    int: number of occurences of sub in string, includes overlaps
    """
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count += 1
        else:
            return count


def import_seeds(seed_file):
    """Import and organize seed information"""

    # read in seed file and organize seed information
    seeds = pd.read_csv(seed_file, sep='\t', header=None).fillna('')
    seeds.columns = ['miRNA family', 'Seed', 'Species with miRNA']
    seeds = seeds.set_index('Seed')
    # seed_window_dict = get_seed_window_dict(list(seeds.index))

    # make a dictionary listing all the species that have a particular microRNA
    seed_to_species = {}
    for row in seeds.iterrows():
        if row[0] not in seed_to_species.keys():
            if row[1]['Species with miRNA'] != '':
                seed_to_species[row[0]] = \
                    row[1]['Species with miRNA'].split(';')
            else:
                seed_to_species[row[0]] = []

    return seed_to_species, seeds


def import_utrs(utr_file):
    """Import and organize utr information"""

    # read in utr file
    UTRS = pd.read_csv(utr_file, sep='\t', header=None).astype(str)
    UTRS.columns = ['Gene', 'Species', 'UTR sequence']

    # sanitize utr sequences
    UTRS['UTR sequence'] = ['X{}X'.format(x.upper().replace('T', 'U')
                                          .replace('\n', ''))
                            for x in UTRS['UTR sequence']]
    UTRS = UTRS.set_index('Gene')

    # make a separate dataframe for the utr sequences of the reference species
    UTRS_REF = UTRS[UTRS['Species'] == config.REF_SPECIES]
    UTRS_REF['UTR sequence'] = [x.replace('-', '')
                                for x in UTRS_REF['UTR sequence']]

    return UTRS, UTRS_REF


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


def is_subtype(type1, type2):
    """Check if site type 1 is contained in site type 2"""
    if type2 == '6mer':
        return True
    elif type1 == '8mer-1a':
        return True
    else:
        return type1 == type2


def get_matches(sp, ref_starts, ref_ends, ref_types,
                this_starts, this_ends, this_types):
    """Check if any of the site locations match those in the reference"""
    l1, l2 = 0, 0
    return_list = [None] * len(ref_starts)
    while (l1 < len(ref_starts)) & (l2 < len(this_starts)):
        start1, end1 = ref_starts[l1], ref_ends[l1]
        start2, end2 = this_starts[l2], this_ends[l2]
        # xxxxxxxx
        #     xxxxxxxx
        if (end1-start2) < config.REQUIRED_OVERLAP:
            l1 += 1
        #     xxxxxxxx
        # xxxxxxxx
        elif (end2-start1) < config.REQUIRED_OVERLAP:
            l2 += 1
        else:
            if is_subtype(this_types[l2], ref_types[l1]):
                return_list[l1] = sp
            l1 += 1
            l2 += 1
    return return_list


def get_site_info(utr_no_gaps, seed):
    """Find site locations and site types for a utr and a miRNA"""

    # (m8, 1a) matches linked to site type
    site_type_dict = {(0, 0): '6mer',
                      (0, 1): '7mer-1a',
                      (1, 0): '7mer-m8',
                      (1, 1): '8mer-1a'}

    # take the reverse complement
    rc_seed = rev_comp(seed)

    # find all the seed match locations
    locs = [m.start() for m
            in re.finditer('(?={})'.format(rc_seed[1:]), utr_no_gaps)]

    # if we find no seed matches, return empty tuples
    if len(locs) == 0:
        return (), (), ()

    # compare the 1a and m8 positions to determine site types
    site_types = [site_type_dict[((utr_no_gaps[x-1] == rc_seed[0]),
                                  (utr_no_gaps[x+6] == 'A'))]
                  for x in locs]

    # find site starts and ends
    site_starts = [x - 1 - int(y in ['8mer-1a', '7mer-m8'])
                   for (x, y) in zip(locs, site_types)]

    site_ends = [x + 5 + int(y in ['8mer-1a', '7mer-1a'])
                 for (x, y) in zip(locs, site_types)]

    return tuple(site_starts), tuple(site_ends), tuple(site_types)


def get_gapped_ends(utr, regex, types):
    """Finds the site ends in the utrs with gaps"""
    type_lengths = {'8mer-1a': 8, '7mer-m8': 7, '7mer-1a': 7, '6mer': 6}
    starts = [m.start() for m in re.finditer('(?={})'.format(regex), utr)]
    ends = []
    for i, s in enumerate(starts):
        if types[i] in ['8mer-1a','7mer-m8']:
            s -= 1
            while utr[s] == '-':
                s -= 1
            starts[i] = s
        length = type_lengths[types[i]] - 1
        end = s+1
        while length > 0:
            if utr[end] != '-':
                length -= 1
            end += 1
        ends.append(end)

    return starts, ends


def find_aligning_species(utr_df, seed, species_with_mirna,
                          num_sites, site_types):
    """
    Finds species with sites in the same location as the reference

    Parameters
    ----------
    utr_df: pandas DataFrame, dataframe of aligned sequences
    seed: string
    species_with_mirna: list of species that contain this miRNA
    num_sites: int, number of sites the reference utr has
    site_types: list of the site types of this miRNA in this gene (in order)

    Output
    ------
    list: list of lists, for each site in the reference, return a list of
        species that have a site in the same location
    """

    # check that there are other species and that this mirna exists in the ref
    if (len(species_with_mirna) == 0) | \
            (config.REF_SPECIES not in species_with_mirna):
        return [''] * num_sites

    utr_df = utr_df.groupby('Species').first()
    ref_utr = utr_df.loc[config.REF_SPECIES]['UTR sequence']

    site = rev_comp(seed)
    regex = '-*'.join(site[1:])

    # find the site locations in the reference, including gaps
    ref_starts, ref_ends = get_gapped_ends(ref_utr, regex, site_types)

    species_list = []

    # loop through species that have this miRNA
    for sp in species_with_mirna:

        # species is not in alignment
        if sp not in utr_df.index:
            continue

        # species is the reference species, don't need to recalculate sites
        # we know they all align to the reference
        elif sp == config.REF_SPECIES:
            species_list.append([sp] * len(ref_starts))

        # look for sites in this species that align to sites in the ref species
        else:
            this_utr = utr_df.loc[sp]['UTR sequence']

            this_types = get_site_info(this_utr.replace('-', ''), seed)[2]

            this_starts, this_ends = get_gapped_ends(this_utr, regex, this_types)

            if len(this_starts) != 0:
                species_list.append(get_matches(sp, ref_starts, ref_ends, site_types,
                                                this_starts, this_ends, this_types))

    # if we didn't find any aligning species, return empty strings
    if len(species_list) == 0:
        return [''] * num_sites

    species_list = zip(*species_list)
    species_list = [[x for x in sp if x is not None] for sp in species_list]

    return species_list


def get_branch_length_score_generic(species_list):
    """Calculate BLS from the generic tree"""

    # get rid of any -1, which indicate the utr does not have the site
    species_list = [x for x in species_list if x != -1]

    # check if any species have the site
    if len(species_list) == 0:
        return 0

    # retreive tree path information from the generic tree
    all_paths = []
    for sp in species_list:
        all_paths += config.SPECIES_TO_PATH[sp]
    all_paths = list(set(list(all_paths)))

    # calculate branch length score
    bls = sum([x[1] for x in all_paths])

    return bls


def get_branch_length_score_specific(species_list, bin):
    """Calculate BLS from a bin-specific tree"""

    # get rid of any -1, which indicate the utr does not have the site
    species_list = [x for x in species_list if x != -1]

    # check if any species have the site
    if len(species_list) == 0:
        return 0

    # retreive tree path information from the bin-specific tree
    all_paths = []
    for sp in species_list:
        all_paths += config.TREES[bin][sp]
    all_paths = list(set(list(all_paths)))

    # calculate branch length score
    bls = sum([x[1] for x in all_paths])

    return bls


def calculate_pct(aligning_species, bin, site_type, seed):
    """From the list of aligning species, calculate PCT"""

    # check if the site type is a 6mer or doesn't exist in the species
    if site_type == '6mer':
        return 0.0, 0.0, 0

    if seed not in config.PARAMS[site_type]:
        return 0.0, 0.0, 0

    # use the helper function to calculate the bls
    bls = get_branch_length_score_specific(aligning_species, bin)

    # retrieve constants
    b0, b1, b2, b3 = config.PARAMS[site_type][seed]

    # calculate PCT as described in the paper
    score = max(0.0, b0 + (b1 / (1.0 + (np.exp(((0.0 - b2) * bls) + b3)))))

    # check if the PCT meets the conservation threshold
    conserved = int(score >= config.CONSERVATION_CUTOFFS[site_type])

    return [bls, score, conserved]
