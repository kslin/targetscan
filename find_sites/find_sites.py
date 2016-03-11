import re
import sys
import time

import concurrent.futures
import numpy as np
import pandas as pd

import config
import find_sites_helpers as fsh
import tasks
import utils

T0 = time.time()

# disable chained assignment warning
pd.options.mode.chained_assignment = None

SEED_FILE = sys.argv[1]
BIN_FILE = sys.argv[2]
UTR_FILE = sys.argv[3]
OUT_FILE = sys.argv[4]


# import seed information
print "Importing seed, bin, and utr files..."
t0 = time.time()

SEED_WINDOW_DICT, SEED_TO_SPECIES, SEEDS = fsh.import_seeds(SEED_FILE)
num_seeds = len(SEEDS)

# Import gene to bin information
BINS = pd.read_csv(BIN_FILE, sep='\t', header=None)
BINS.columns = ['Gene', 'BLS', 'Bin']
BINS = BINS.set_index('Gene')

# import UTR information
UTRS, UTRS_REF = fsh.import_utrs(UTR_FILE)

print "{} seconds".format(time.time() - t0)


# Find sites
print 'Finding sites...'
t0 = time.time()

# create a dataframe that lists every miRNA/utr combination possible
gene_list = list(UTRS_REF.index) * num_seeds
sequence_list = list(UTRS_REF['UTR sequence']) * num_seeds
bin_list = list(BINS.loc[list(UTRS_REF.index)]['Bin'].fillna(1)) * num_seeds
num_genes = len(UTRS_REF)
seed_list = []
family_list = []
for row in SEEDS.iterrows():
    seed_list += [row[0]] * num_genes
    family_list += [row[1]['miRNA family']] * num_genes

site_df = pd.DataFrame({
    'Gene ID': gene_list,
    'UTR sequence': sequence_list,
    'Bin':  bin_list,
    'Seed': seed_list,
    'miRNA family': family_list
    })

# filter out rows where the miRNA seed does not appear in the UTR sequence
site_df['Num sites'] = [utils.occurrences(utr,
                                          utils.reverse_complement(seed[:-1]))
                        for (utr, seed)
                        in zip(site_df['UTR sequence'], site_df['Seed'])]
site_df = site_df[site_df['Num sites'] > 0]

# group the data by gene and pass to a helper function to get site information
groups = site_df.groupby('Gene ID')
num_genes = len(groups)
site_info = []

# if indicated by the user, use the paralle implementation
if config.FUTURES:
    executor = concurrent.futures.ProcessPoolExecutor()
    futures = []
    for gene, group in groups:
        futures.append(executor.submit(tasks.get_all_site_info, gene,
                                       group, UTRS.loc[gene], SEED_TO_SPECIES,
                                       SEED_WINDOW_DICT))

        # add a sleep to prevent the executor from getting clogged
        time.sleep(0.0001)

    # as jobs are completed, add the information
    for future in concurrent.futures.as_completed(futures):
        site_info += future.result()

    executor.shutdown()

# otherwise, use the nonparallel implementation
else:
    for gene, group in groups:
        site_info += tasks.get_all_site_info(gene, group, UTRS.loc[gene],
                                             SEED_TO_SPECIES, SEED_WINDOW_DICT)

# convert site information into a dataframe and add column names
site_info = pd.DataFrame(site_info)
site_info.columns = ['Gene ID', 'miRNA family', 'Seed', 'Site start',
                     'Site end', 'Site type', 'PCT', 'Conserved',
                     'Branch length score', 'Aligning species']

site_info['UTR sequence'] = [x.replace('X', '')
                             for x in UTRS_REF.loc[site_info['Gene ID']]
                             ['UTR sequence']]
site_info.loc[:, 'UTR BLS'] = list(BINS.loc[site_info['Gene ID']]['BLS'])

new_sites, new_starts = list(site_info['Site type']), list(site_info['Site start'])
for i in range(len(new_sites)):
    if new_starts[i] == config.TOO_CLOSE - 1:
        site = new_sites[i]
        if site == '8mer-1a':
            new_sites[i] = '7mer-1a'
            new_starts[i] += 1
        elif site == '7mer-m8':
            new_sites[i] = '6mer'
            new_starts[i] += 1

site_info['Site type'] = new_sites
site_info['Site start'] = new_starts
site_info = site_info[site_info['Site start'] >= config.TOO_CLOSE]
print "{} seconds".format(time.time() - t0)


# Write to file
print "Writing outfile..."
t0 = time.time()

site_info.to_csv(OUT_FILE, sep='\t', index=False)

print "{} seconds".format(time.time() - t0)

print "Total time elapsed: {}".format(time.time() - T0)
