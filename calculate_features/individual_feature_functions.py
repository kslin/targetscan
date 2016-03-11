import os
import re
import shlex
import subprocess
import sys

import numpy as np
import pandas as pd

import config
import utils


def get_rnaplfold_data(gene,utr):

    # sanitize name of file so we don't break the shell
    gene_name = shlex.split(gene)[0]

    cwd = os.getcwd()
    os.chdir(config.RNAPLFOLD_FOLDER)

    mytempfile = 'temp_{}.fa'.format(gene_name)
    with open(mytempfile,'wb') as f:
        f.write('>{}\n{}'.format(gene_name,utr))


    mycall = 'RNAplfold -L {} -W {} -u 20 < {}'.format(min(40,len(utr)),min(80,len(utr)),mytempfile)
    subprocess.call([mycall],shell=True,stdout=subprocess.PIPE)
    lunp_file = '{}_lunp'.format(gene_name)

    rnaplfold_data = pd.read_csv(lunp_file,sep='\t',header=1).set_index(' #i$')

    os.remove(mytempfile)
    os.remove(lunp_file)
    os.remove('{}_dp.ps'.format(gene_name))

    os.chdir(cwd)

    return rnaplfold_data


def calculate_threep_score(mirna,utr,site_type,site_start):
    mirna_3p = mirna[8:]
    utr_5p = utils.reverse_complement(utr[max(0,site_start-15):site_start - int(site_type in ['6mer','7mer-1a'])])
    scores = np.empty((len(utr_5p)+1,len(mirna_3p)+1))
    scores.fill(np.nan)
    possible_scores = [0]
    for i,nt1 in enumerate(utr_5p):
        for j,nt2 in enumerate(mirna_3p):
            if nt1 == nt2:
                new_score = 0.5 + 0.5*((j>3)&(j<8))
                if np.isnan(scores[i,j]) == False:
                    new_score += scores[i,j]
                    scores[i+1,j+1] = new_score
                    possible_scores.append(new_score)
                else:
                    scores[i+1,j+1] = new_score - max(0,(abs(i-j)-2)*0.5)
            else:
                scores[i+1,j+1] = float('NaN')
    return np.nanmax(possible_scores)

def calculate_local_au(utr,site_type,site_start,site_end):
    up_site_adder = int(site_type not in ['7mer-m8','8mer-1a'])
    upstream = utr[max(0,site_start-30):site_start]
    upstream = [int(x in ['A','U']) for x in upstream]
    upweights = [1.0/(x+1+up_site_adder) for x in range(len(upstream))][::-1]

    down_site_adder = int(site_type in ['7mer-1a','8mer-1a'])
    downstream = utr[site_end:min(len(utr),site_end+30)]
    downstream = [int(x in ['A','U']) for x in downstream]
    downweights = [1.0/(x+1+down_site_adder) for x in range(len(downstream))]

    return (np.dot(upstream,upweights) + np.dot(downstream,downweights))/(sum(upweights)+sum(downweights))

def calculate_min_dist(site_start,site_end,airs_subdf):
    airs_subdf = airs_subdf[airs_subdf['AIR end'] >= site_end]
    min_dists = [min(site_start,x-site_end) for x in airs_subdf['AIR end'].values]
    weighted = np.dot(min_dists,airs_subdf['AIR ratio'].values)/float(np.sum(airs_subdf['AIR ratio'].values))
    if weighted <= 1:
            return 0
    return np.log10(weighted)

def calculate_weighted_utr_length(site_end,airs_subdf):
    airs_subdf = airs_subdf[airs_subdf['AIR end'] >= site_end]
    weighted = np.dot(airs_subdf['AIR end'],airs_subdf['AIR ratio'])/float(np.sum(airs_subdf['AIR ratio']))
    if weighted <= 1:
            return 0
    return np.log10(weighted)


def calculate_weighted_num_off6mers(off6m_locs,site_end,airs_subdf):
    if len(off6m_locs) == 0:
        return 0
    else:
        airs_subdf = airs_subdf[airs_subdf['AIR end'] >= site_end]
        airs_subdf.loc[:,'counts'] = [len([loc for loc in off6m_locs if loc < end]) for end in airs_subdf['AIR end']]
        weighted = np.dot(airs_subdf['counts'],airs_subdf['AIR ratio'])/float(np.sum(airs_subdf['AIR ratio']))
        if weighted <= 1:
            return 0
        return np.log10(weighted)

# utr = 'AUUGCACUUUAUAUGUCCAGCCUCUUCCUCAGUCCCCCAAACCAAAGCUACACAGCCAGAUUUCUCAAGCAGUCUCAACUCCAGUCCCUCAUCUCACCCUUACUAUUGCUCUUGCUUUCCAGUUUGCUUUUGAUUUGCAUCUUCUCACUAGUAAAACUGCCUUCCCUUUGUUCCUUAUUUUCUGUUUUUUCUCUAGAGAGGUACAGUUGUAAGUCAGAGUUAAUAUAAUAGGGCCUGUGAAAACAGAGGCUUUUGCAUUGUCUCUUGACAUCAGAAGUUACAAUAGGCAUAUGGGCAAAAUGGUGUAGCAGGCUCACUGGCCGUUUGUUUUUUAAACACAUUUUCACAAGUUUUUGAGACACUGGAUUUCUUUAAUUAAAAAAAAAAUGCCAAGAAACAUUAUUUAUACAGGGUUGAUUGCUUUCAUGUUGUUAUUCUGUACCCUAUAGUAGCCUCCAUGAGAAUCUGGUAUUUCUUGCUGCUUGGAACUACUUUGCAGUGAUUACUUGGUUGCAGUCCAAGUACUCUCGUUUAGUCUGAGCCUGGAGAUGUUCUAGACUUGCUUCUCCCACCUCUGAGAUUAGGACAGGAAAAAUGUGAAAUUUCCCAAUUACAGGAUUAUACGGUACCAUCACAUCAUUUGUGGAAAUUGGGGUGACUGUAUAGCUGGGAUUGGGCUAAGGACUGUGGUCUUAUCUGUCCACAUACAGCCAAAAUGCCUAUCCAGAAAUCCAGUUCGUUGGAAAGGAAAAUUGGUACUCCUGUGCCACAGGGGUUCCAGAAAAGGGAAGUCACUUUACCUUGCGGUGGUGGGAUCCUGAUGUCUUUCAUCCAUUUGUAGUAAAAGCUGGUAAAGCUUUUCUUACUCCUGGUUCCCUACCAGUAUUUCUAAACAUGUCGCACUUUCUCCACAGgcaugugguuuugaccuuuuuuucaaucuucuagaaagggaacggaagcagaagugggacaucgagggcucugcuguccucugcgcuggguguggaaugcugcugcaccugucccuucugcuggcucagggaagugucuucuugcccacauuucuguggggaaagguuuuuaauccucUGAUGCUUCCAUCUUCCUGUUUAGGCCAUGUGCCCAGAAACCUGGACUGAUCUUUCUUUAAUAGUGAACCCCUGGGCCACUGAAGAGUAACAUGGCUCCACUGGACACAAAAGAGGGAUGGAAUCAACAGGCAGGGGGCCUUUUAUAAGCCUUAGGAAAAGAAAAUGAAACUAUUUCAUCUUUGGACUUUUCAAUACUAUUGGAGUGAUUUUUUUCUUUCUAAACAGGGAAAAUAAUGUUACAAAAGCAUCUUUUUUGUUAUUUGUUUGCAUCCCUCCCCCACACCCUGGUGUUUUAAAAUGAAGAAAAAAAACCAUCACCUUUUGUACAAAAACUCUUAAUGAUUAAAAAACAAACAAAACAUA'
# print len(utr)
# utr = utr.upper()
# # seed = 'AAAGUGC'
# mirna = 'UAAAGUGCUGACAGUGCAGAU'

# print calculate_threep_score(mirna, utr,'7mer-1a',792)


# def calculate_SA(site_type,site_start,lunp_file):
#     site_start_for_SA = site_start
#     if site_type in ['6mer','7mer-1a']:
#         site_start_for_SA -= 1
#     with open(lunp_file,'rb') as f:
#         for i,line in enumerate(f):
#             if i == site_start_for_SA+9:
#                 assert (int(line.split('\t')[0]) == site_start_for_SA+8), "{},{}".format(line.split('\t')[0],site_start_for_SA+8)
#                 sa_score = np.log10(float(line.split('\t')[14]))
#                 return sa_score
#     return 0

