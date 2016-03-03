import time

import numpy as np
import pandas as pd

import config
import bin_helper_functs


def get_median_bls(i,gene,group):
    # if i %1000 == 0:
    #     print i
    group = group.set_index('Species ID')
    if config.REF_SPECIES not in group.index:
        return [gene,0,0]

    ref_utr = group.loc[config.REF_SPECIES][['UTR sequence']].values[0]
    group = group.drop(config.REF_SPECIES)

    if len(group) == 0:
        return [gene,0,0]

    ref_utr = ref_utr.replace('-','$')
    num_gap = ref_utr.count('$')

    utr_list = zip(*[
        bin_helper_functs.compare_seqs(x, ref_utr, group.index[i]) for i, x in enumerate(group['UTR sequence'])
    ])

    subdf = pd.DataFrame({
        'nts': utr_list,
        'count': [1] * len(ref_utr)
    })
    subdf = subdf.groupby('nts').agg({'count':np.nansum})
    subdf['bls'] = [bin_helper_functs.get_branch_length_score(x) for x in subdf.index]
    subdf = subdf.sort_values(by=['bls'])

    bls = bin_helper_functs.get_median_binned_list(subdf['count'], subdf['bls'], len(ref_utr)-num_gap, num_gap)

    return [gene, bls, bin_helper_functs.get_bin(bls)]
