import os

import bin_helpers


REF_SPECIES = os.environ.get('REF_SPECIES', '9606')
print("ref_species: {}".format(REF_SPECIES))

BL_THRESHOLDS = [0, 1.21207417, 2.17396073, 2.80215414, 3.26272822,
                 3.65499277, 4.01461968, 4.40729032, 4.90457274, 5.78196252]

TREE_FILE_GENERIC = '../PCT_parameters/Tree.generic.txt'
SPECIES_TO_PATH = bin_helpers.parse_tree(TREE_FILE_GENERIC, REF_SPECIES)
