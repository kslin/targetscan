import re

import pandas as pd

import config
import individual_feature_functions
import utils


AIRS = pd.read_csv(config.AIRS_FILE,sep='\t',header=None)
AIRS.columns = ['Gene ID','AIR start','AIR end','AIR ratio']
AIRS = AIRS.set_index('Gene ID')

def calculate_all_features(i, inputs):
    gene, mirna, utr, seed, site_type, site_start, site_end = inputs
    rc_seed = utils.reverse_complement(seed)
    off6 = rc_seed[:-1]
    off6m_locs = list(set([m.start() for m in re.finditer(off6, utr)]) - set([m.start() for m in re.finditer(rc_seed, utr)]))
    airs_subdf = AIRS.loc[[gene]]

    three_p_score = individual_feature_functions.calculate_threep_score(mirna, utr, site_type, site_start)
    local_au_score = individual_feature_functions.calculate_local_au(utr, site_type, site_start, site_end)
    min_dist_score = individual_feature_functions.calculate_min_dist(site_start, site_end, airs_subdf)
    utr_length_score = individual_feature_functions.calculate_weighted_utr_length(site_end, airs_subdf)
    off6mer_score = individual_feature_functions.calculate_weighted_num_off6mers(off6m_locs, site_end, airs_subdf)

    return [i, three_p_score, local_au_score, min_dist_score, utr_length_score, off6mer_score]