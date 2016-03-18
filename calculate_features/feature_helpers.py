import os
import shlex
from string import maketrans
import subprocess
import sys

import numpy as np
import pandas as pd

import config


def rev_comp(seq):
    """
    Parameters:
    ==========
    seq: string, sequence to get reverse complement of

    Returns:
    =======
    float: reverse complement of seq in all caps
    """
    seq = seq.upper()
    intab = "AUCG"
    outtab = "UAGC"
    trantab = maketrans(intab, outtab)
    seq = seq[::-1]
    seq = seq.translate(trantab)
    return seq


def get_rnaplfold_data(gene, utr):
    """
    Run RNAplfold and get pairing probabilities for a utr

    Parameters
    ----------
    gene: string, name of gene

    utr: string, utr sequence

    Output
    ------
    pandas DataFrame: pairing probabilities at each position
    """

    # sanitize name of file so we don't break the shell
    gene_name = shlex.split(gene)[0]

    # navigate to the folder for RNAplfold data
    cwd = os.getcwd()
    os.chdir(config.RNAPLFOLD_FOLDER)

    # write sequence to a temporary file
    mytempfile = 'temp_{}.fa'.format(gene_name)
    with open(mytempfile, 'wb') as f:
        f.write('>{}\n{}'.format(gene_name, utr))

    # call RNAplfold
    length = min(40, len(utr))
    window = min(80, len(utr))
    mycall = 'RNAplfold -L {} -W {} -u 20 < {}'.format(length, window,
                                                       mytempfile)
    subprocess.call([mycall], shell=True, stdout=subprocess.PIPE)
    lunp_file = '{}_lunp'.format(gene_name)

    # read data and convert to a dataframe
    rnaplfold_data = pd.read_csv(lunp_file, sep='\t',
                                 header=1).set_index(' #i$')

    os.remove(mytempfile)
    os.remove(lunp_file)
    os.remove('{}_dp.ps'.format(gene_name))

    os.chdir(cwd)

    return rnaplfold_data


def calculate_threep_score(mirna, utr, site_type, site_start):
    """
    Calculate the three-prime pairing score

    Parameters
    ----------
    mirna: string, miRNA sequence

    utr: string, utr sequence

    site_type: string, site type

    site_start: int, start of site

    Output
    ------
    float: 3' pairing score
    """
    # get the 3' region of the mirna and the corresponding utr seq
    mirna_3p = mirna[8:]
    trailing = utr[max(0, site_start-15):
                   site_start - int(site_type in ['6mer', '7mer-1a'])]
    utr_5p = rev_comp(trailing)

    # initiate array for dynamic programming search
    scores = np.empty((len(utr_5p) + 1, len(mirna_3p) + 1))
    scores.fill(np.nan)
    possible_scores = [0]

    # fill in array
    for i, nt1 in enumerate(utr_5p):
        for j, nt2 in enumerate(mirna_3p):
            if nt1 == nt2:
                new_score = 0.5 + 0.5 * ((j > 3) & (j < 8))
                if np.isnan(scores[i, j]) == False:
                    new_score += scores[i, j]
                    scores[i + 1, j + 1] = new_score
                    possible_scores.append(new_score)
                else:
                    offset_penalty = max(0, (abs(i - j) - 2) * 0.5)
                    scores[i + 1, j + 1] = new_score - offset_penalty
            else:
                scores[i + 1, j + 1] = float('NaN')

    return np.nanmax(possible_scores)


def calculate_local_au(utr, site_type, site_start, site_end):
    """
    Calculate the local AU score

    Parameters
    ----------
    utr: string, utr sequence

    site_type: string, site type

    site_start: int, start of site

    site_end: int, end of site

    Output
    ------
    float: local AU score
    """
    # find A, U and weights upstream of site
    up_site_adder = int(site_type not in ['7mer-m8', '8mer-1a'])
    upstream = utr[max(0, site_start - 30): site_start]
    upstream = [int(x in ['A', 'U']) for x in upstream]
    upweights = [1.0 / (x + 1 + up_site_adder)
                 for x in range(len(upstream))][::-1]

    # find A,U and weights downstream of site
    down_site_adder = int(site_type in ['7mer-1a', '8mer-1a'])
    downstream = utr[site_end:min(len(utr), site_end + 30)]
    downstream = [int(x in ['A', 'U']) for x in downstream]
    downweights = [1.0 / (x + 1 + down_site_adder)
                   for x in range(len(downstream))]

    weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
    total = float(sum(upweights) + sum(downweights))

    return weighted / total


def calculate_min_dist(site_start, site_end, airs_subdf):
    """
    Calculate the min dist score

    Parameters
    ----------
    site_start: int, start of site

    site_end: int, end of site

    airs_subdf: pandas DataFrame, alternative isoform ratios for this gene

    Output
    ------
    float: min dist score
    """
    # only get alternative isoforms that contain this site
    airs_subdf = airs_subdf[airs_subdf['AIR end'] >= site_end]
    min_dists = [min(site_start, x - site_end)
                 for x in airs_subdf['AIR end'].values]
    total = float(np.sum(airs_subdf['AIR ratio'].values))
    weighted = np.dot(min_dists, airs_subdf['AIR ratio'].values) / total

    # check if we can take the log
    if weighted <= 1:
        return 0

    return np.log10(weighted)


def calculate_weighted_utr_length(site_end, airs_subdf):
    """
    Calculate the utr length score

    Parameters
    ----------
    site_end: int, end of site

    airs_subdf: pandas DataFrame, alternative isoform ratios for this gene

    Output
    ------
    float: utr length score
    """
    # only get alternative isoforms that contain this site
    airs_subdf = airs_subdf[airs_subdf['AIR end'] >= site_end]
    total = float(np.sum(airs_subdf['AIR ratio']))
    weighted = np.dot(airs_subdf['AIR end'], airs_subdf['AIR ratio']) / total

    # check we can take the log
    if weighted <= 1:
        return 0

    return np.log10(weighted)


def calculate_weighted_num_off6mers(off6m_locs, site_end, airs_subdf):
    """
    Calculate the weighted number of offset 6mers score

    Parameters
    ----------
    off6m_locs: list of ints, locations of offset 6mers

    site_end: int, end of site

    airs_subdf: pandas DataFrame, alternative isoform ratios for this gene

    Output
    ------
    float: offset 6mer score
    """
    if len(off6m_locs) == 0:
        return 0
    else:
        airs_subdf = airs_subdf[airs_subdf['AIR end'] >= site_end]
        airs_subdf.loc[:, 'counts'] = [len([loc for loc in off6m_locs
                                            if loc < end])
                                       for end in airs_subdf['AIR end']]
        total = float(np.sum(airs_subdf['AIR ratio']))
        weighted = np.dot(airs_subdf['counts'], airs_subdf['AIR ratio']) / total

        # check if we can take the log
        if weighted <= 1:
            return 0

        return np.log10(weighted)
