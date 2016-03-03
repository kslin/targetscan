import os

import pandas as pd

import find_sites_helpers

FUTURES = os.environ.get('FUTURES', False)
REF_SPECIES = '9606'
TOO_CLOSE = 15

TREE_FILE_GENERIC = '../PCT_parameters/Tree.generic.txt'
SPECIES_TO_PATH = find_sites_helpers.parse_tree(TREE_FILE_GENERIC, REF_SPECIES)

TREE_PATH = '../PCT_parameters/'
TREES = {}
for i in range(10):
	TREES[i+1] = find_sites_helpers.parse_tree(TREE_PATH + 'Tree.bin_{:02}.txt'.format(i+1), REF_SPECIES)

PARAMS = {}
for site_type in ['8mer-1a','7mer-m8','7mer-1a']:
	params = {}
	with open(TREE_PATH + site_type.replace('-','_')+'_PCT_parameters.txt','rb') as f:
		f.next()
		for line in f:
			line = line.split('\t')
			params[line[0]] = (float(line[1]),float(line[2]),float(line[3]),float(line[4]))

	PARAMS[site_type] = params

CONSERVATION_CUTOFFS = {'8mer-1a': 1.8, '7mer-m8': 2.8, '7mer-1a': 3.6}