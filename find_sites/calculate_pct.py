import math

import numpy as np
import pandas as pd

import config


def get_branch_length_score2(species_list,bin):
    species_list = [x for x in species_list if x != -1]
    if len(species_list) == 0:
        return 0
    all_paths = []
    for sp in species_list:
        all_paths += config.TREES[bin][sp]
    all_paths = list(set(list(all_paths)))
    bls = sum([x[1] for x in all_paths])
    return bls

def calculate_pct(aligning_species,bin,site_type,seed):
    if site_type == '6mer':
        return 0.0

    if seed not in config.PARAMS[site_type]:
        return 0.0

    bls = get_branch_length_score2(aligning_species,bin)
    b0,b1,b2,b3 = config.PARAMS[site_type][seed]

    score = max(0.0, b0 + (b1/(1.0 + (np.exp(((0.0-b2)*bls) + b3)))))

    conserved = int(score >= config.CONSERVATION_CUTOFFS[site_type])

    # return [score, conserved]
    return score

mir_data = pd.read_csv('../../infiles/miR_Family_info_sample.txt',sep='\t',header=None)
mir_data.columns = ['miRNA family','seed','species']
mir_data = mir_data.set_index('miRNA family')

ts_output = pd.read_csv('../../outfiles/targetscan_sites_tiny_output.txt',sep='\t').fillna('')
ts_output = ts_output[ts_output['species_ID'] == 9606]

ts_pct = pd.read_csv('../../outfiles/targetscan_tiny_output.txt',sep='\t').set_index('Gene ID')
ts_pct = ts_pct[ts_pct['Species ID'] == 9606]
ts_pct = ts_pct[['miRNA family','Site Type','UTR start','PCT contribution']]

bins = pd.read_csv('../../outfiles/bins_tiny_outfile.txt',sep='\t',header=None)
bins.columns = ['Gene ID','BLS','bin']
bins = bins.set_index('Gene ID')


# print ts_pct[ts_pct['PCT contribution'] > 0]



# ts_output = ts_output[ts_output['a_Gene_ID'] == 'CDC2L6']

# print ts_output


# ts_output = ts_output[ts_output['Species_in_this_group_with_this_site_type'] != '']
# print ts_output

for row in ts_output.iterrows():
    # print row
    if len(row[1]['Species_in_this_group']) > 0:
        family = row[1]['miRNA_family_ID']
        seed = mir_data.loc[family]['seed']
        start = row[1]['UTR_start']
        subdf = ts_pct[ts_pct['miRNA family'] == family]
        subdf = subdf[subdf['UTR start'] == start].loc[[row[1]['a_Gene_ID']]]
        pct = subdf.iloc[0]['PCT contribution']
        bin = bins.loc[row[1]['a_Gene_ID']]['bin']
        aligning = row[1]['Species_in_this_group'].split(' ')
        # print row[1]['a_Gene_ID'],row[1]['miRNA_family_ID'],row[1]['UTR_start'],aligning, bin, row[1]['Site_type'] 
        print pct - calculate_pct(aligning,bin,row[1]['Site_type'],seed)


