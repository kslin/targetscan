import re
import sys
import time

import numpy as np
import pandas as pd

import config
import utils


# disable chained assignment warning
pd.options.mode.chained_assignment = None

MIRNA_FILE = sys.argv[1]
UTR_FILE = sys.argv[2]
OUT_FILE = sys.argv[3]

## Import miRNA and seed information ##
t0 = time.time()
print "Processing miRNAs..."

mirnas = pd.read_csv(MIRNA_FILE,sep='\t',header=None).fillna('')
mirnas.columns = ['miRNA family','Seed','Mirbase ID','miRNA sequence','Species_with_miRNA']
seed_to_species = {}
for row in mirnas.iterrows():
    if row[1]['Seed'] not in seed_to_species.keys():
        if row[1]['Species_with_miRNA'] != '':
            seed_to_species[row[1]['Seed']] = row[1]['Species_with_miRNA'].split(';')
        else:
            seed_to_species[row[1]['Seed']] = []

assert len(mirnas) == len(list(set(list(mirnas['miRNA sequence'])))), "Error: duplicate mirna sequences in input file"

SEEDS = list(set(list(mirnas['Seed'])))
SEED_WINDOW_DICT = find_sites_helpers.get_seed_window_dict(SEEDS)
num_mirnas = len(mirnas)

print "{} seconds".format(time.time()-t0)

## Import UTR information ##
t0 = time.time()
print "Processing UTRs..."
genes = []
sequences = []
species = []
with open(UTR_FILE,'rb') as utr_file:
    for line in utr_file:
        line = line.split('\t')
        if int(line[1]) in ALL_SPECIES:
            genes.append(line[0])
            species.append(int(line[1]))
            sequences.append('X' + line[2].replace('-','').upper().replace('T','U') + 'X')

UTRS = pd.DataFrame({'Gene':genes,'Species':species,'UTR sequence': sequences})
UTRS = UTRS[UTRS['Species'] == config.REF_SPECIES].drop('Species',1)
UTRS = UTRS.set_index('Gene')

print "{} seconds".format(time.time()-t0)


## Find UTRs with seed matches ##
t0 = time.time()
print "Extracting miRNA and UTR information..."
gene_list = list(UTRS.index)*num_mirnas
sequence_list = list(UTRS['UTR sequence'])*num_mirnas
num_genes = len(UTRS)
seed_list = []
mirseq_list = []
family_list = []
mirbase_list = []
for row in mirnas.iterrows():
    seed = row[1]['Seed']
    seed_list += [seed] * num_genes
    mirseq_list += [row[1]['miRNA sequence']] * num_genes
    family_list += [row[1]['miRNA family']] * num_genes
    mirbase_list += [row[1]['Mirbase ID']] * num_genes

site_df = pd.DataFrame({
    'Gene ID':gene_list,
    'Seed':seed_list,
    'miRNA sequence':mirseq_list,
    'miRNA family':family_list,
    'UTR sequence':sequence_list,
    'Mirbase ID':mirbase_list
    })

site_df.loc[:,'num sites'] = [utils.occurrences(utr,utils.reverse_complement(seed[:-1])) for (utr,seed) in zip(site_df['UTR sequence'],site_df['Seed'])]
site_df = site_df[site_df['num sites'] > 0]
site_df['index'] = range(len(site_df))
site_df = site_df.set_index('index')

print "{} seconds".format(time.time()-t0)

## Find locations and site types of all sites ##
t0 = time.time()
print "Finding sites..."

merge_df = pd.DataFrame([find_sites_helpers.get_site_info(utr,seed,SEED_WINDOW_DICT[seed]) for utr,seed in zip(site_df['UTR sequence'],site_df['Seed'])],index = site_df.index)
merge_df.columns = ['Site start','Site end','Site type']

site_df = pd.concat([site_df,merge_df],axis=1,join_axes=[site_df.index])

site_df['UTR sequence'] = [x.replace('X','') for x in site_df['UTR sequence']]

print "{} seconds".format(time.time()-t0)

## Write to file ##
t0 = time.time()
print "Writing outfile..."

site_df.to_csv(OUT_FILE,sep='\t',index=False)

print "{} seconds".format(time.time()-t0)

print "Total time elapsed: {}".format(time.time()-T0)

