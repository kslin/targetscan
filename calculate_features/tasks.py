import re
import os
import shlex
import shutil
import subprocess
import sys

import numpy as np
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
    utr = group.iloc[0]['UTR sequence']
    utr_bls = group.iloc[0]['UTR BLS']
    rnaplfold_data = individual_feature_functions.get_rnaplfold_data(gene,utr)
    data = []
    for row in group.iterrows():
        seed, family, site_type, site_start, site_end, pct, conserved, bls = row[1]['Seed'], row[1]['miRNA family'], row[1]['Site type'], row[1]['Site start'], row[1]['Site end'], row[1]['PCT'], row[1]['Conserved'], row[1]['Branch length score']
        mirna_subdf = mirna_df.loc[[family]]
        mirnas = list(mirna_subdf['miRNA sequence'])
        mirbases = list(mirna_subdf['Mirbase ID'])
        num_mirnas = len(mirnas)
        rc_seed = utils.reverse_complement(seed)
        off6 = rc_seed[:-1]
        off6m_locs = list(set([m.start() for m in re.finditer(off6, utr)]) - set([m.start() for m in re.finditer(rc_seed, utr)]))

        # depends on rest of miRNA
        three_p_scores = [individual_feature_functions.calculate_threep_score(mirna, utr, site_type, site_start) for mirna in mirnas]
        mirna1As, mirna1Cs, mirna1Gs = zip(*[[int(mirna[0] == nt) for nt in ['A','C','G']] for mirna in mirnas])

        # only depends on seed and site
        local_au_score = individual_feature_functions.calculate_local_au(utr, site_type, site_start, site_end)
        min_dist_score = individual_feature_functions.calculate_min_dist(site_start, site_end, airs_subdf)
        utr_length_score = individual_feature_functions.calculate_weighted_utr_length(site_end, airs_subdf)
        off6mer_score = individual_feature_functions.calculate_weighted_num_off6mers(off6m_locs, site_end, airs_subdf)
        mirna8A, mirna8C, mirna8G = (int(seed[-1] == nt) for nt in ['A','C','G'])
        site8A, site8C, site8G = (int(utr[site_start-1] == nt) if site_type in ['6mer','7mer-1a'] else 0 for nt in ['A','C','G'])

        site_start_for_SA = site_start - int(site_type in ['6mer','7mer-1a'])
        if (site_start_for_SA + 8) not in rnaplfold_data.index:
            sa_score = 0
        else:
            raw_sa_score = rnaplfold_data.loc[site_start_for_SA+8]['14']
            if raw_sa_score <= 0:
                sa_score = 0
            else:
                sa_score = np.log10(raw_sa_score)

        combined = [[gene,family,mirbase,mirna,seed,site_type,site_start,site_end,three_p_score,local_au_score,min_dist_score,utr_length_score,off6mer_score,sa_score,mirna1A,mirna1C,mirna1G,mirna8A,mirna8C,mirna8G,site8A,site8C,site8G,pct,conserved,bls,utr_bls] for three_p_score,mirna,mirbase,mirna1A,mirna1C,mirna1G in zip(three_p_scores,mirnas,mirbases,mirna1As,mirna1Cs,mirna1Gs)]

        data += combined

    return data
