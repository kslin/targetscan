import ast
import re
from string import maketrans
# from subprocess import call
import sys
import time

import concurrent.futures
import numpy as np
import pandas as pd

import config
from tasks import calculate_all_features
import utils


T0 = time.time()

# disable chained assignment warning
pd.options.mode.chained_assignment = None

assert len(sys.argv) == 5, "See README for usage"
ORF_file, TARGET_file, BIN_FILE, OUT_file = sys.argv[1:]

## IMPORT ALL FILES AND SAVE TO DATAFRAMES ##
t0 = time.time()
print "Importing files..."

ORFS = pd.read_csv(ORF_file,sep='\t',header=None).astype(str)
ORFS.columns = ['Gene ID','Species ID','ORF sequence']
ORFS = ORFS[ORFS['Species ID'] == config.REF_SPECIES]
ORFS['ORF sequence'] = [x.replace('-','').upper().replace('T','U') for x in ORFS['ORF sequence']]
ORFS = ORFS.set_index('Gene ID')

print '{} seconds'.format(time.time()-t0)

# # BINS = pd.read_csv(BIN_FILE,sep='\t')

TARGETS = pd.read_csv(TARGET_file,sep='\t').astype(str)
TARGETS = TARGETS[TARGETS['Species ID'] == config.REF_SPECIES].drop('Species ID',1)
TARGETS[['Site start','Site end']] = TARGETS[['Site start','Site end']].astype(int)

## Add ORF data ##
t0 = time.time()
print "Adding ORF data..."

TARGETS.loc[:,'ORF sequence'] = list(ORFS.loc[TARGETS['Gene ID']]['ORF sequence'])
TARGETS.loc[:,'ORF length'] = [np.log10(len(orf)) for orf in TARGETS['ORF sequence']]
zipped = zip(list(TARGETS['ORF sequence']),list(TARGETS['Seed']))
TARGETS.loc[:,'ORF 8mers'] = [orf.count(utils.reverse_complement(seed) + 'A') for (orf,seed) in zipped]
TARGETS = TARGETS.drop('ORF sequence',1)

print '{} seconds'.format(time.time()-t0)

## Add TA and SPS data ##
t0 = time.time()
print "Adding TA and SPS data..."

TA_SPS = pd.read_csv(config.TA_SPS_FILE,sep='\t')
TA_SPS = TA_SPS.set_index('Seed region')
TA_SPS = TA_SPS.loc[TARGETS['Seed']]
TA_SPS.loc[:,'Site type'] = list(TARGETS['Site type'])
zipped = zip(list(TA_SPS['SPS (8mer and 7mer-m8)']),list(TA_SPS['SPS (7mer-1a and 6mer)']),list(TA_SPS['Site type']))
TA_SPS.loc[:,'SPS'] = [a if stype in ['8mer-1a','7mer-m8'] else b for (a,b,stype) in zipped]

TARGETS.loc[:,'TA'] = list(TA_SPS['TA'])
TARGETS.loc[:,'SPS'] = list(TA_SPS['SPS'])

TARGETS['index'] = range(len(TARGETS))
TARGETS = TARGETS.set_index('index')

print '{} seconds'.format(time.time()-t0)

## Calculate other features ##
t0 = time.time()
print "Calculating features..."

zipped = zip(list(TARGETS['Gene ID']),list(TARGETS['miRNA sequence']),list(TARGETS['UTR sequence']),
    list(TARGETS['Seed']),list(TARGETS['Site type']),list(TARGETS['Site start']),list(TARGETS['Site end']))
TARGETS = TARGETS.drop(['UTR sequence'],1)

data = np.zeros((len(TARGETS),6))

# run parallel version
if config.FUTURES:
    executor = concurrent.futures.ProcessPoolExecutor()
    futures = []
    for i,inputs in enumerate(zipped):
    	futures.append(executor.submit(calculate_all_features, i, inputs))
    	time.sleep(0.00001)

    print time.time()-T0

    i = 0
    for future in concurrent.futures.as_completed(futures):
    	data[i,:] = future.result()
    	i += 1

    executor.shutdown()

# otherwise, run non-parallel version
else:
    for i,inputs in enumerate(zipped):
    	data[i,:] = calculate_all_features(i, inputs)

merge_df = pd.DataFrame(data)
merge_df.columns = ['index', 'Threep score', 'Local AU score', 'Min dist', 'UTR length score', 'Off6m score']
merge_df = merge_df.set_index('index')
TARGETS = pd.concat([TARGETS, merge_df], axis=1, join='inner')

print '{} seconds'.format(time.time()-t0)

## Write to file ##
t0 = time.time()
print "Writing to file..."

TARGETS.to_csv(OUT_file,sep='\t',index=False)

print '{} seconds'.format(time.time()-t0)

print 'Total time elapsed: {}'.format(time.time()-T0)



