import os


FUTURES = os.environ.get('FUTURES', False)

BL_THRESHOLDS = [0,1.21207417,2.17396073,2.80215414,3.26272822,3.65499277,4.01461968,4.40729032,4.90457274,5.78196252]
REF_SPECIES = '9606'
RNAPLFOLD_FOLDER = './RNAplfold_temp'
TA_SPS_FILE = 'TA_SPS_by_seed_region.txt'

## Tiny data ##

# AIRS_FILE = '../../infiles/All_cell_lines.AIRs_small.txt'
# MIRNA_FILE = '../../infiles/miR_for_context_scores.sample.txt'
# TARGET_FILE = '../../outfiles/sites_tiny_outfile.txt' 
# ORF_FILE = '../../infiles/ORF_Sequences_sample.txt'
# OUT_FILE = '../../outfiles/tiny_output.txt'

## Small data ##

AIRS_FILE = '../../infiles/All_cell_lines.AIRs.txt'
# MIRNA_FILE = '../../infiles/kathy_mirna_file.txt'
# TARGET_FILE = '../../outfiles/sites_small_outfile.txt' 
# ORF_FILE = '../../infiles/ORF_Sequences_human.txt'
# OUT_FILE = '../../outfiles/small_output.txt'

## Big data ##

# AIRS_FILE = '../../infiles/All_cell_lines.AIRs.txt'
# MIRNA_FILE = '../../infiles/kathy_mirna_file.txt'
# TARGET_FILE = '../../outfiles/sites_big_outfile.txt' 
# ORF_FILE = '../../infiles/ORF_Sequences_human.txt'
# OUT_FILE = '../../outfiles/big_output.txt'