import sys
import time

import numpy as np
import pandas as pd

import config
import bls_helpers


def get_bls(gene, group):
    """
    Given a group of aligned sequences, find the bls at every position.

    Parameters
    ----------
    gene: string, name of gene

    group: pandas DataFrame, aligned sequences indexed by species ID

    Output
    ------
    list of bls values at each position in the utr
    """

    group = group.set_index('Species ID')

    # if the reference species is not in the aligned sequences, return empty list
    if config.REF_SPECIES not in group.index:
        return []

    # otherwise, isolate the reference utr sequence and drop from the rest
    ref_utr = group.loc[config.REF_SPECIES][['UTR sequence']].values[0]
    group = group.drop(config.REF_SPECIES)

    # if the reference utr was the only one in the alignment, return empty list
    if len(group) == 0:
        return []


    # for each position in the ref sequence, find aligning species
    utr_list = zip(*[
        bls_helpers.compare_seqs(x, ref_utr, group.index[i])
        for i, x in enumerate(group['UTR sequence'])
    ])

    # group by aligning species so we only have to calculate each bls once
    patterns = list(set(utr_list))
    bls_dict = {x: bls_helpers.get_branch_length_score(x)
                for x in patterns}

    # construct dataframe with values
    subdf = pd.DataFrame({
        'nts': list(ref_utr),
        'pattern': utr_list
    })

    # get rid of rows with gaps
    subdf = subdf[subdf['nts'] != '-']

    # add in bls
    subdf['bls'] = [bls_dict[x] for x in subdf['pattern']]

    return [gene] + list(subdf['bls'])
