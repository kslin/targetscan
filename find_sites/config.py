import os

import pandas as pd

import find_sites_helpers as fsh


# define global variables
FUTURES = os.environ.get('FUTURES', False)
REF_SPECIES = '9606'
TOO_CLOSE = 14
CONSERVATION_CUTOFFS = {'8mer-1a': 1.8, '7mer-m8': 2.8, '7mer-1a': 3.6}

# specify paths to parameter files
TREE_FILE_GENERIC = '../PCT_parameters/Tree.generic.txt'
TREE_PATH = '../PCT_parameters/'

# parse the generic tree
SPECIES_TO_PATH = fsh.parse_tree(TREE_FILE_GENERIC, REF_SPECIES)

# parse bin-specific trees
TREES = {}
for i in range(10):
    TREES[i+1] = fsh.parse_tree(TREE_PATH + 'Tree.bin_{:02}.txt'.format(i+1),
                                REF_SPECIES)

# parse constants for calculating PCT
PARAMS = {}
for site_type in ['8mer-1a', '7mer-m8', '7mer-1a']:
    params = {}
    tree_file = TREE_PATH + site_type.replace('-', '_')+'_PCT_parameters.txt'
    with open(tree_file, 'rb') as f:
        f.next()
        for line in f:
            line = line.split('\t')
            params[line[0]] = tuple([float(line[i]) for i in [1, 2, 3, 4]])

    PARAMS[site_type] = params
