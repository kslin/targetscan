import re
import sys
import time

import concurrent.futures
import numpy as np
import pandas as pd

import config
import find_sites_helpers
import tasks
import utils

T0 = time.time()

# disable chained assignment warning
pd.options.mode.chained_assignment = None

SEED_FILE = sys.argv[1]
BIN_FILE = sys.argv[2]
UTR_FILE = sys.argv[3]
OUT_FILE = sys.argv[4]


## Import seed information ##
t0 = time.time()
print "Processing seed file..."

SEEDS = pd.read_csv(SEED_FILE,sep='\t',header=None).fillna('')
SEEDS.columns = ['miRNA family','Seed','Species with miRNA']
SEEDS = SEEDS.set_index('Seed')
num_seeds = len(SEEDS)
SEED_WINDOW_DICT = find_sites_helpers.get_seed_window_dict(list(SEEDS.index))

SEED_TO_SPECIES = {}
for row in SEEDS.iterrows():
    if row[0] not in SEED_TO_SPECIES.keys():
        if row[1]['Species with miRNA'] != '':
            SEED_TO_SPECIES[row[0]] = row[1]['Species with miRNA'].split(';')
        else:
            SEED_TO_SPECIES[row[0]] = []

print "{} seconds".format(time.time()-t0)

## Import UTR information ##
t0 = time.time()
print "Processing UTRs..."

UTRS = pd.read_csv(UTR_FILE,sep='\t',header=None).astype(str)
UTRS.columns = ['Gene','Species','UTR sequence']
UTRS['UTR sequence'] = ['X' + x.upper().replace('T','U').replace('\n','') + 'X' for x in UTRS['UTR sequence']]
UTRS = UTRS.set_index('Gene')

UTRS_REF = UTRS[UTRS['Species'] == config.REF_SPECIES]
UTRS_REF['UTR sequence'] = [x.replace('-','') for x in UTRS_REF['UTR sequence']]

print "{} seconds".format(time.time()-t0)

## Import gene to bin information ##
t0 = time.time()
print 'Importing bin information...'

BINS = pd.read_csv(BIN_FILE,sep='\t',header=None)
BINS.columns = ['Gene','BLS','Bin']
BINS = BINS.set_index('Gene')


gene_list = list(UTRS_REF.index)*num_seeds
sequence_list = list(UTRS_REF['UTR sequence'])*num_seeds
# bin_list = list(BINS.loc[list(UTRS_REF.index)]['Bin'].fillna(1))*num_seeds
num_genes = len(UTRS_REF)
seed_list = []
family_list = []
for row in SEEDS.iterrows():
    seed_list += [row[0]]*num_genes
    family_list += [row[1]['miRNA family']]*num_genes

site_df = pd.DataFrame({
    'Gene ID':gene_list,
    'UTR sequence':sequence_list,
    # 'Bin': bin_list,
    'Seed':seed_list,
    'miRNA family':family_list
    })

site_df.loc[:,'Num sites'] = [utils.occurrences(utr,utils.reverse_complement(seed[:-1])) for (utr,seed) in zip(site_df['UTR sequence'],site_df['Seed'])]
site_df = site_df[site_df['Num sites'] > 0]

# zipped = zip(site_df['Gene ID'],site_df['Seed'],site_df['Num sites'],site_df['UTR sequence'],site_df['miRNA family'])
groups = site_df.groupby('Gene ID')
num_genes = len(groups)
site_info = []
if config.FUTURES:
    executor = concurrent.futures.ProcessPoolExecutor()
    futures = []
    for i,(gene,group) in enumerate(groups):
        futures.append(executor.submit(tasks.get_all_site_info,gene,group,UTRS.loc[gene],SEED_TO_SPECIES,SEED_WINDOW_DICT))
    # for i,(gene,seedm8,num_sites,utr_no_gaps,family) in enumerate(zipped):
    #     futures.append(executor.submit(tasks.get_all_site_info,gene,family,seedm8,UTRS.loc[gene],SEED_TO_SPECIES[seedm8],num_sites,utr_no_gaps,SEED_WINDOW_DICT[seedm8]))
        if (i % 1000) == 0:
            print '{}/{}'.format(i,num_genes)

    # for future in tqdm(concurrent.futures.as_completed(futures), total=len(UTRS), unit='UTR', leave=True):
    for future in concurrent.futures.as_completed(futures):
        site_info += future.result()
        # time.sleep(0.00001)

    print "{} seconds".format(time.time()-t0)

    executor.shutdown()
else:
    # for gene, utr_group in tqdm(UTRS, unit='UTR', leave=True):
    for i,(gene,group) in enumerate(groups):
    # for i,(gene,seedm8,num_sites,utr_no_gaps,family) in enumerate(zipped):
        # site_info += tasks.get_all_site_info(gene,family,seedm8,UTRS.loc[gene],SEED_TO_SPECIES[seedm8],num_sites,utr_no_gaps,SEED_WINDOW_DICT[seedm8])
        site_info += tasks.get_all_site_info(gene,group,UTRS.loc[gene],SEED_TO_SPECIES,SEED_WINDOW_DICT)
        if (i%1000) == 0:
            print '{}/{}'.format(i,num_genes)
    print "{} seconds".format(time.time()-t0)

# site_df.loc[:,'Aligning species'] = aligning_species
site_info = pd.DataFrame(site_info)
site_info.columns = ['Gene ID','miRNA family','Seed','Site start','Site end','Site type','Branch length score']
site_info.loc[:,'UTR sequence'] = [x.replace('X','') for x in UTRS_REF.loc[site_info['Gene ID']]['UTR sequence']]
site_info.loc[:,'UTR BLS'] = list(BINS.loc[site_info['Gene ID']]['BLS'])

## Write to file ##
t0 = time.time()
print "Writing outfile..."

site_info.to_csv(OUT_FILE,sep='\t',index=False)

print "{} seconds".format(time.time()-t0)

print "Total time elapsed: {}".format(time.time()-T0)


# print site_df.head()
sys.exit()
