import os

import bls_helpers


FUTURES = os.environ.get('FUTURES', False)

BL_THRESHOLDS = [0, 1.21207417, 2.17396073, 2.80215414, 3.26272822,
                 3.65499277, 4.01461968, 4.40729032, 4.90457274, 5.78196252]
REF_SPECIES = '9606'
TREE_FILE_GENERIC = '../PCT_parameters/Tree.generic.txt'
SPECIES_TO_PATH = bls_helpers.parse_tree(TREE_FILE_GENERIC, REF_SPECIES)
