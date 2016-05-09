import os

import pandas as pd


FUTURES = os.environ.get('FUTURES', False)

BL_THRESHOLDS = [0, 1.21207417, 2.17396073, 2.80215414,
                 3.26272822, 3.65499277, 4.01461968,
                 4.40729032, 4.90457274, 5.78196252]
REF_SPECIES = '9606'
RNAPLFOLD_FOLDER = './RNAplfold_temp'
TA_SPS_FILE = 'TA_SPS_by_seed_region.txt'

# AIRS_FILE = '../../infiles/All_cell_lines.AIRs_small.txt'
AIRS_FILE = '../../infiles/All_cell_lines.AIRs.txt'

# read in AIRS file
AIRS = pd.read_csv(AIRS_FILE, sep='\t', header=None)
AIRS.columns = ['Gene ID', 'AIR start', 'AIR end', 'AIR ratio']
AIRS = AIRS.set_index('Gene ID')
