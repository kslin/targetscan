import sys
import time

import numpy as np
import pandas as pd

import config
import bin_helpers


def get_median_bls(gene, group):
    """
    Given a group of aligned sequences, find the median bls and bin.

    Parameters
    ----------
    gene: string, name of gene

    group: pandas DataFrame, aligned sequences indexed by species ID

    Output
    ------
    [string, float, int]: [gene, branch length score, bin]
    """

    group = group.set_index('Species ID')

    # if the reference species is not in the aligned sequences, return 0
    if config.REF_SPECIES not in group.index:
        return [gene, 0, 0]

    # otherwise, isolate the reference utr sequence and drop from the rest
    ref_utr = group.loc[config.REF_SPECIES][['UTR sequence']].values[0]
    group = group.drop(config.REF_SPECIES)

    # if the reference utr was the only one in the alignment, return 0
    if len(group) == 0:
        return [gene, 0, 0]

    # replace gaps with another character and count the gaps
    ref_utr = ref_utr.replace('-', '$').upper()
    num_gap = ref_utr.count('$')

    # for each position in the ref sequence, find aligning species
    utr_list = zip(*[
        bin_helpers.compare_seqs(x, ref_utr, group.index[i])
        for i, x in enumerate(group['UTR sequence'])
    ])

    # group by aligning species so we only have to calculate each bls once
    subdf = pd.DataFrame({
        'nts': utr_list,
        'count': [1] * len(ref_utr)
    })

    subdf = subdf.groupby('nts').agg({'count': np.nansum})
    subdf['bls'] = [bin_helpers.get_branch_length_score(x)
                    for x in subdf.index]
    subdf = subdf.sort_values(by=['bls'])

    # use the helper function to get the branch length score
    bls = bin_helpers.get_median_binned_list(subdf['count'],
                                             subdf['bls'],
                                             len(ref_utr)-num_gap, num_gap)

    bls = round(bls, 8)

    return [gene, bls, bin_helpers.get_bin(bls)]
