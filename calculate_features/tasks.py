import re
import os
import shlex
import subprocess
import sys

import numpy as np
import pandas as pd

import feature_helpers as feat


def calculate_all_features(gene, group, mirna_df, rnaplfold_temp, airs_subdf=None):
    """
    Given a gene and a group of miRNAs with targets in that gene,
        calculate all the features

    Parameters
    ----------
    gene: string, name of gene

    group: pandas DataFrame, site information

    mirna_df: pandas DataFrame, miRNA information

    Output
    ------
    list of lists: list of features for each miRNA
    """

    # get gene-specific data
    utr = group.iloc[0]['UTR sequence']
    utr_bls = group.iloc[0]['UTR BLS']

    # run RNAplfold
    rnaplfold_data = feat.get_rnaplfold_data(gene, utr, rnaplfold_temp)

    data = []
    for row in group.iterrows():
        seed, site_type, site_start, site_end = (row[1]['Seed'],
                                                 row[1]['Site type'],
                                                 row[1]['Site start'],
                                                 row[1]['Site end'])

        family, pct, conserved, bls = (row[1]['miRNA family'],
                                       row[1]['PCT'],
                                       row[1]['Conserved'],
                                       row[1]['Branch length score'])

        # get list of different miRNAs for a given seed
        mirna_subdf = mirna_df.loc[[family]]
        mirnas = list(mirna_subdf['miRNA sequence'])
        mirbases = list(mirna_subdf['Mirbase ID'])
        num_mirnas = len(mirnas)

        # calculate offset 6mer counts
        rc_seed = feat.rev_comp(seed)
        off6 = rc_seed[:-1]

        # find where offset 6mer locations are
        off6m_locs_total = [m.start() for m in re.finditer(off6, utr)]

        # find other site locations, so we can eliminate cases where the
        # offset 6mer is a subsite of a larger site
        off6mer_locs_overlaps = [m.start() for m in re.finditer(rc_seed, utr)]

        # find the difference in these two sets
        off6m_locs = list(set(off6m_locs_total) - set(off6mer_locs_overlaps))

        # features that depend on rest of miRNA
        three_p_scores = [feat.calculate_threep_score(mirna, utr,
                                                      site_type, site_start)
                          for mirna in mirnas]
        mirna1As, mirna1Cs, mirna1Gs = zip(*[[int(mirna[0] == nt)
                                              for nt in ['A', 'C', 'G']]
                                             for mirna in mirnas])

        # only depends on seed and site
        local_au_score = feat.calculate_local_au(utr, site_type,
                                                 site_start, site_end)

        # if no airs_subdf given, calculate these features for whole UTR
        if airs_subdf is None:
            utr_length = len(utr)
            min_dist_score = np.log10(min(site_start, utr_length - site_end + 1))
            utr_length_score = np.log10(utr_length)
            off6mer_score = np.log10(len(off6m_locs) + 1)
        else:
            min_dist_score = feat.calculate_min_dist_withAIRs(site_start,
                                                     site_end,
                                                     airs_subdf)

            utr_length_score = feat.calculate_weighted_utr_length_withAIRs(site_end,
                                                                  airs_subdf)

            off6mer_score = feat.calculate_weighted_num_off6mers_withAIRs(off6m_locs,
                                                                 site_end,
                                                                 airs_subdf)

        mirna8A, mirna8C, mirna8G = tuple([int(seed[-1] == nt)
                                           for nt in ['A', 'C', 'G']])

        site8A, site8C, site8G = tuple([int(utr[site_start-1] == nt)
                                        if site_type in ['6mer', '7mer-1a']
                                        else 0 for nt in ['A', 'C', 'G']])

        # use the rnaplfold data to calculate the site accessibility
        site_start_for_SA = site_start - int(site_type in ['6mer', '7mer-1a'])
        if (site_start_for_SA + 8) not in rnaplfold_data.index:
            sa_score = 0
        else:
            raw_sa_score = rnaplfold_data.loc[site_start_for_SA + 8]['14']
            if raw_sa_score <= 0:
                sa_score = 0
            else:
                sa_score = np.log10(raw_sa_score)

        zipped = zip(three_p_scores, mirnas, mirbases,
                     mirna1As, mirna1Cs, mirna1Gs)

        combined = [[gene, family, mirbase, mirna, seed,
                     site_type, site_start, site_end,
                     three_p_score, local_au_score, min_dist_score,
                     utr_length_score, off6mer_score, sa_score,
                     mirna1A, mirna1C, mirna1G, mirna8A, mirna8C, mirna8G,
                     site8A, site8C, site8G, pct, conserved, bls, utr_bls]
                    for (three_p_score, mirna, mirbase,
                         mirna1A, mirna1C, mirna1G)
                    in zipped]

        data += combined

    return data
