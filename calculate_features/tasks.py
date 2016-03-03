import re

import pandas as pd

import config
import individual_feature_functions
import utils


AIRS = pd.read_csv(config.AIRS_FILE,sep='\t',header=None)
AIRS.columns = ['Gene ID','AIR start','AIR end','AIR ratio']
AIRS = AIRS.set_index('Gene ID')

# def calculate_all_features(i, inputs):
#     gene, mirna, utr, seed, site_type, site_start, site_end = inputs
#     rc_seed = utils.reverse_complement(seed)
#     off6 = rc_seed[:-1]
#     off6m_locs = list(set([m.start() for m in re.finditer(off6, utr)]) - set([m.start() for m in re.finditer(rc_seed, utr)]))
#     airs_subdf = AIRS.loc[[gene]]

#     three_p_score = individual_feature_functions.calculate_threep_score(mirna, utr, site_type, site_start)
#     local_au_score = individual_feature_functions.calculate_local_au(utr, site_type, site_start, site_end)
#     min_dist_score = individual_feature_functions.calculate_min_dist(site_start, site_end, airs_subdf)
#     utr_length_score = individual_feature_functions.calculate_weighted_utr_length(site_end, airs_subdf)
#     off6mer_score = individual_feature_functions.calculate_weighted_num_off6mers(off6m_locs, site_end, airs_subdf)

#     return [i, three_p_score, local_au_score, min_dist_score, utr_length_score, off6mer_score]

def calculate_all_features(gene, group, mirna_df):
    airs_subdf = AIRS.loc[[gene]]
    data = []
    for row in group.iterrows():
        seed, family, utr, site_type, site_start, site_end, bls, utr_bls = row[1]['Seed'], row[1]['miRNA family'], row[1]['UTR sequence'], row[1]['Site type'], row[1]['Site start'], row[1]['Site end'], row[1]['Branch length score'], row[1]['UTR BLS']
        mirna_subdf = mirna_df.loc[[family]]
        mirnas = list(mirna_subdf['miRNA sequence'])
        mirbases = list(mirna_subdf['Mirbase ID'])
        num_mirnas = len(mirnas)
        rc_seed = utils.reverse_complement(seed)
        off6 = rc_seed[:-1]
        off6m_locs = list(set([m.start() for m in re.finditer(off6, utr)]) - set([m.start() for m in re.finditer(rc_seed, utr)]))

        three_p_scores = [individual_feature_functions.calculate_threep_score(mirna, utr, site_type, site_start) for mirna in mirnas]

        local_au_score = individual_feature_functions.calculate_local_au(utr, site_type, site_start, site_end)
        min_dist_score = individual_feature_functions.calculate_min_dist(site_start, site_end, airs_subdf)
        utr_length_score = individual_feature_functions.calculate_weighted_utr_length(site_end, airs_subdf)
        off6mer_score = individual_feature_functions.calculate_weighted_num_off6mers(off6m_locs, site_end, airs_subdf)

        combined = [[gene,family,mirbase,mirna,seed,site_type,site_start,site_end,three_p_score,local_au_score,min_dist_score,utr_length_score,off6mer_score,bls,utr_bls] for three_p_score,mirna,mirbase in zip(three_p_scores,mirnas,mirbases)]

        data += combined

    return data
