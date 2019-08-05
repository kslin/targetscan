#!/usr/bin/python2

from optparse import OptionParser
import os
import sys
import time

import concurrent.futures
import numpy as np
import pandas as pd

import feature_helpers
from tasks import calculate_all_features


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--mirna_file", dest="MIRNA_FILE", help="file with miRNA sequences")
    parser.add_option("--site_file", dest="SITE_FILE", help="output from find_sites")
    parser.add_option("--ta_sps_file", dest="TA_SPS_FILE", help="TA and SPS by seed region")
    parser.add_option("--orf_file", dest="ORF_FILE", help="ORF sequences in tab delimited format")
    parser.add_option("--rnaplfold_temp", dest="RNAPLFOLD_TEMP", help="temp folder to write RNAplfold outputs")
    parser.add_option("--out", dest="OUT_FILE", help="where to write bin output")
    parser.add_option("--ref_species", dest="REF_SPECIES", help="reference species", default='9606')
    parser.add_option("--airs_file", dest="AIRS_FILE", default=None, help="affected isoform ratios")
    parser.add_option("--futures", dest="FUTURES", help="if true, run in parallel", default=False, action='store_true')

    (options, args) = parser.parse_args()


    T0 = time.time()

    # disable chained assignment warning
    pd.options.mode.chained_assignment = None

    # import target and mirna data
    t0 = time.time()
    print "Adding target and miRNA data..."

    # read in miRNA information
    MIRNAS = pd.read_csv(options.MIRNA_FILE, sep='\t', header=None).astype(str)
    MIRNAS.columns = ['miRNA family',
                      'Species ID',
                      'Mirbase ID',
                      'miRNA sequence']
    MIRNAS = MIRNAS[MIRNAS['Species ID'] == options.REF_SPECIES]
    MIRNAS = MIRNAS.set_index('miRNA family')

    # read in target information and set variable types
    TARGETS = pd.read_csv(options.SITE_FILE, sep='\t').astype(str)

    TARGETS[['Site start', 'Site end']] = \
        TARGETS[['Site start', 'Site end']].astype(int)

    TARGETS[['UTR BLS', 'Branch length score']] = \
        TARGETS[['UTR BLS', 'Branch length score']].astype(float)

    # read in AIRs table if given
    if options.AIRS_FILE is not None:
        AIRS = pd.read_csv(options.AIRS_FILE, sep='\t', header=None)
        AIRS.columns = ['Gene ID', 'AIR start', 'AIR end', 'AIR ratio']
        AIRS = AIRS.set_index('Gene ID')
    else:
        AIRS = None

    print '{} seconds\n'.format(time.time() - t0)

    # calculate features
    t0 = time.time()
    print "Calculating features..."

    # group data by gene
    groups = TARGETS.groupby('Gene ID')
    num_genes = len(groups)

    # make a folder for RNAPLFOLD data
    if os.path.isdir(options.RNAPLFOLD_TEMP) is False:
        os.mkdir(options.RNAPLFOLD_TEMP)

    data = []

    # run calculations in parallel if indicated
    if options.FUTURES:
        print("Running parallel version")
        executor = concurrent.futures.ProcessPoolExecutor()
        futures = []
        for i, (gene, group) in enumerate(groups):

            # get the isoform ratio table for this gene if supplied
            if AIRS is not None:
                airs_subdf = AIRS.loc[[gene]]
            else:
                airs_subdf = None

            # calculate features in parallel
            futures.append(executor.submit(
                calculate_all_features,
                gene,
                group,
                MIRNAS,
                options.RNAPLFOLD_TEMP,
                airs_subdf
            ))

            # add sleep so we don't overwhelm the executor
            time.sleep(0.0001)

        # aggregate data as it is computed
        for future in concurrent.futures.as_completed(futures):
            data += future.result()

        executor.shutdown()

    # otherwise, run non-parallel version
    else:
        for i, (gene, group) in enumerate(groups):

            # get the isoform ratio table for this gene if supplied
            if options.AIRS_FILE is not None:
                airs_subdf = AIRS.loc[[gene]]
            else:
                airs_subdf = None

            # calculate features
            data += calculate_all_features(gene, group, MIRNAS, options.RNAPLFOLD_TEMP, airs_subdf)
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

    print '{} seconds\n'.format(time.time() - t0)

    # add TA and SPS data
    t0 = time.time()
    print "Adding TA and SPS data..."

    TA_SPS = pd.read_csv(options.TA_SPS_FILE, sep='\t')
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

    print '{} seconds\n'.format(time.time() - t0)

    # add ORF length and ORF 8mer count
    t0 = time.time()
    print "Adding ORF data..."

    ORFS = pd.read_csv(options.ORF_FILE, sep='\t', header=None).astype(str)
    ORFS.columns = ['Gene ID', 'Species ID', 'ORF sequence']
    ORFS = ORFS[ORFS['Species ID'] == options.REF_SPECIES]
    ORFS['ORF sequence'] = [x.replace('-', '').upper().replace('T', 'U')
                            for x in ORFS['ORF sequence']]
    ORFS = ORFS.set_index('Gene ID')

    orf_list = list(ORFS.loc[TARGETS['Gene ID']].fillna('')['ORF sequence'])
    TARGETS['ORF length'] = [np.log10(len(orf)) if len(orf) > 0 else 0
                             for orf in orf_list]
    zipped = zip(orf_list, list(TARGETS['Seed']))
    TARGETS['ORF 8mers'] = [orf.count(feature_helpers.rev_comp(seed) + 'A')
                            for (orf, seed) in zipped]

    print '{} seconds\n'.format(time.time() - t0)

    # write to file
    t0 = time.time()
    print "Writing to file..."

    TARGETS.to_csv(options.OUT_FILE, sep='\t', index=False, float_format='%.6f')

    print '{} seconds\n'.format(time.time() - t0)

    print 'Total time for calculating features: {}\n'.format(time.time() - T0)
