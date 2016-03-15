import ast
import os
import re
from string import maketrans
import subprocess
import sys
import time

import concurrent.futures
import numpy as np
import pandas as pd

import config
from tasks import calculate_all_features
import utils


def main(MIRNA_FILE, TARGET_FILE, ORF_FILE, OUT_FILE):
    T0 = time.time()

    # disable chained assignment warning
    pd.options.mode.chained_assignment = None

    # import target and mirna data
    t0 = time.time()
    print "Adding target and miRNA data..."

    # read in miRNA information
    MIRNAS = pd.read_csv(MIRNA_FILE, sep='\t', header=None).astype(str)
    MIRNAS.columns = ['miRNA family',
                      'Species ID',
                      'Mirbase ID',
                      'miRNA sequence']
    MIRNAS = MIRNAS[MIRNAS['Species ID'] == config.REF_SPECIES]
    MIRNAS = MIRNAS.set_index('miRNA family')

    # read in target information and set variable types
    TARGETS = pd.read_csv(TARGET_FILE, sep='\t').astype(str)

    TARGETS[['Site start', 'Site end']] = \
        TARGETS[['Site start', 'Site end']].astype(int)

    TARGETS[['UTR BLS', 'Branch length score']] = \
        TARGETS[['UTR BLS', 'Branch length score']].astype(float)

    print '{} seconds'.format(time.time() - t0)

    # calculate features
    t0 = time.time()
    print "Calculating features..."

    # group data by gene
    groups = TARGETS.groupby('Gene ID')
    num_genes = len(groups)

    # make a folder for RNAPLFOLD data
    if os.path.isdir(config.RNAPLFOLD_FOLDER) == False:
        os.mkdir(config.RNAPLFOLD_FOLDER)

    data = []

    # run calculations in parallel if indicated
    if config.FUTURES:
        executor = concurrent.futures.ProcessPoolExecutor()
        futures = []
        for i, (gene, group) in enumerate(groups):
            futures.append(executor.submit(calculate_all_features,
                                           gene, group, MIRNAS))

            # add sleep so we don't overwhelm the executor
            time.sleep(0.0001)

        # aggregate data as it is computed
        for future in concurrent.futures.as_completed(futures):
            data += future.result()

        executor.shutdown()

    # otherwise, run non-parallel version
    else:
        for i, (gene, group) in enumerate(groups):
            data += calculate_all_features(gene, group, MIRNAS)
            if (i % 1000) == 0:
                print '{}/{}'.format(i, num_genes)

    TARGETS = pd.DataFrame(data)
    TARGETS.columns = ['Gene ID', 'miRNA family', 'Mirbase ID',
                       'miRNA sequence', 'Seed', 'Site type',
                       'Site start', 'Site end', 'Threep score',
                       'Local AU score', 'Min dist score',
                       'UTR length score', 'Off6m score', 'SA',
                       'siRNA 1A', 'siRNA 1C', 'siRNA 1G',
                       'siRNA 8A', 'siRNA 8C', 'siRNA 8G',
                       'site 8A', 'site 8C', 'site 8G',
                       'PCT', 'Conserved', 'Branch length score', 'UTR BLS']

    os.rmdir(config.RNAPLFOLD_FOLDER)
    print '{} seconds'.format(time.time() - t0)

    # add TA and SPS data
    t0 = time.time()
    print "Adding TA and SPS data..."

    TA_SPS = pd.read_csv(config.TA_SPS_FILE, sep='\t')
    TA_SPS = TA_SPS.set_index('Seed region')
    TA_SPS = TA_SPS.loc[TARGETS['Seed']]
    TA_SPS['Site type'] = list(TARGETS['Site type'])
    zipped = zip(list(TA_SPS['SPS (8mer and 7mer-m8)']),
                 list(TA_SPS['SPS (7mer-1a and 6mer)']),
                 list(TA_SPS['Site type']))
    TA_SPS['SPS'] = [a if stype in ['8mer-1a', '7mer-m8'] else b
                     for (a, b, stype) in zipped]

    TARGETS['TA'] = list(TA_SPS['TA'])
    TARGETS['SPS'] = list(TA_SPS['SPS'])

    print '{} seconds'.format(time.time() - t0)

    # add ORF length and ORF 8mer count
    t0 = time.time()
    print "Adding ORF data..."

    ORFS = pd.read_csv(ORF_FILE, sep='\t', header=None).astype(str)
    ORFS.columns = ['Gene ID', 'Species ID', 'ORF sequence']
    ORFS = ORFS[ORFS['Species ID'] == config.REF_SPECIES]
    ORFS['ORF sequence'] = [x.replace('-', '').upper().replace('T', 'U')
                            for x in ORFS['ORF sequence']]
    ORFS = ORFS.set_index('Gene ID')

    orf_list = list(ORFS.loc[TARGETS['Gene ID']].fillna('')['ORF sequence'])
    TARGETS['ORF length'] = [np.log10(len(orf)) if len(orf) > 0 else 0
                             for orf in orf_list]
    zipped = zip(orf_list, list(TARGETS['Seed']))
    TARGETS['ORF 8mers'] = [orf.count(utils.reverse_complement(seed) + 'A')
                            for (orf, seed) in zipped]

    print '{} seconds'.format(time.time() - t0)

    # write to file
    t0 = time.time()
    print "Writing to file..."

    TARGETS.to_csv(OUT_FILE, sep='\t', index=False)

    print '{} seconds'.format(time.time() - t0)

    print 'Total time elapsed: {}'.format(time.time() - T0)


if __name__ == '__main__':

    if len(sys.argv) != 5:
        print 'See README for usage'
        sys.exit()

    # Get file of aligned UTRS and the file for writing output
    MIRNA_FILE, TARGET_FILE, ORF_FILE, OUT_FILE = sys.argv[1:]

    # run main code
    main(MIRNA_FILE, TARGET_FILE, ORF_FILE, OUT_FILE)
