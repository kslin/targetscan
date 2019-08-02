import os

import pandas as pd


FUTURES = os.environ.get('FUTURES', False)
REF_SPECIES = os.environ.get('REF_SPECIES', '9606')
print("ref_species: {}".format(REF_SPECIES))

RNAPLFOLD_FOLDER = './RNAplfold_temp'
TA_SPS_FILE = 'TA_SPS_by_seed_region.txt'

# AIRS_FILE = '../../infiles/All_cell_lines.AIRs_small.txt'
AIRS_FILE = '../../infiles/All_cell_lines.AIRs.txt'

# read in AIRS file
AIRS = pd.read_csv(AIRS_FILE, sep='\t', header=None)
AIRS.columns = ['Gene ID', 'AIR start', 'AIR end', 'AIR ratio']
AIRS = AIRS.set_index('Gene ID')
