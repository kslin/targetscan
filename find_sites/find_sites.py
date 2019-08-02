#!/usr/bin/python2

from optparse import OptionParser
import os
import sys
import time

import concurrent.futures
import numpy as np
import pandas as pd

import find_sites_helpers as fsh
import tasks


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--seed_file", dest="SEED_FILE", help="file with miRNA seeds")
    parser.add_option("--bin_file", dest="BIN_FILE", help="output from calculate_bins")
    parser.add_option("--utr3_file", dest="UTR3_FILE", help="3\'UTR sequences in tab delimited format")
    parser.add_option("--tree_path", dest="TREE_PATH", help="path to all phylogenetic trees in newick format")
    parser.add_option("--out", dest="OUT_FILE", help="where to write bin output")
    parser.add_option("--ribosome_shadow", dest="RIBOSOME_SHADOW", type=int, default=14, help="length of ribosome shadow")
    parser.add_option("--ref_species", dest="REF_SPECIES", help="reference species", default='9606')
    parser.add_option("--futures", dest="FUTURES", help="if true, run in parallel", default=False, action='store_true')

    (options, args) = parser.parse_args()

    T0 = time.time()

    # disable chained assignment warning
    pd.options.mode.chained_assignment = None

    # parse the generic tree
    SPECIES_TO_PATH = fsh.parse_tree(os.path.join(options.TREE_PATH, 'Tree.generic.txt'), options.REF_SPECIES)

    # parse bin-specific trees
    TREES = {}
    for i in range(10):
        TREES[i+1] = fsh.parse_tree(os.path.join(options.TREE_PATH, 'Tree.bin_{:02}.txt'.format(i+1)),
                                    options.REF_SPECIES)

    # import seed information
    print "Importing seed, bin, and utr files..."
    t0 = time.time()

    SEED_TO_SPECIES, SEEDS = fsh.import_seeds(options.SEED_FILE)
    num_seed = len(SEEDS)

    # parse constants for calculating PCT
    PCT_PARAMS = []
    for site_type in ['8mer-1a', '7mer-m8', '7mer-1a']:
        temp = pd.read_csv(os.path.join(options.TREE_PATH, site_type.replace('-', '_')+'_PCT_parameters.txt'), sep='\t')
        temp['Site_type'] = site_type

        PCT_PARAMS.append(temp)
    PCT_PARAMS = pd.concat(PCT_PARAMS).set_index(['Seed', 'Site_type'])

    # Import gene to bin information
    BINS = pd.read_csv(options.BIN_FILE, sep='\t', header=None)
    BINS.columns = ['Gene', 'BLS', 'Bin']
    BINS = BINS.set_index('Gene')

    # import UTR information
    UTRS, UTRS_REF = fsh.import_utrs(options.UTR3_FILE, options.REF_SPECIES)

    print "{} seconds\n".format(time.time() - t0)

    # Find sites
    print 'Finding sites...'
    t0 = time.time()

    # create a dataframe that lists every miRNA/utr combination possible
    gene_list = list(UTRS_REF.index) * num_seed
    sequence_list = list(UTRS_REF['UTR sequence']) * num_seed
    bin_list = list(BINS.loc[list(UTRS_REF.index)]['Bin'].fillna(1)) * num_seed
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
    site_df['Num sites'] = [fsh.occurrences(utr, fsh.rev_comp(seed[:-1]))
                            for (utr, seed)
                            in zip(site_df['UTR sequence'], site_df['Seed'])]
    site_df = site_df[site_df['Num sites'] > 0]

    # group the data by gene and pass to a helper function to get site info
    groups = site_df.groupby('Gene ID')
    num_genes = len(groups)
    site_info = []

    # if indicated by the user, use the parallel implementation
    if options.FUTURES:
        executor = concurrent.futures.ProcessPoolExecutor()
        futures = []
        for gene, group in groups:
            futures.append(executor.submit(
                tasks.get_all_site_info,
                gene,
                group,
                UTRS.loc[[gene]],
                SEED_TO_SPECIES,
                options.REF_SPECIES,
                SPECIES_TO_PATH,
                TREES,
                PCT_PARAMS
            ))

            # add a sleep to prevent the executor from getting clogged
            time.sleep(0.0001)

        # as jobs are completed, add the information
        for future in concurrent.futures.as_completed(futures):
            site_info += future.result()

        executor.shutdown()

    # otherwise, use the nonparallel implementation
    else:
        for gene, group in groups:
            site_info += tasks.get_all_site_info(
                gene, group,
                UTRS.loc[[gene]],
                SEED_TO_SPECIES,
                options.REF_SPECIES,
                SPECIES_TO_PATH,
                TREES,
                PCT_PARAMS
            )

    # convert site information into a dataframe and add column names
    site_info = pd.DataFrame(site_info)
    site_info.columns = ['Gene ID', 'miRNA family', 'Seed', 'Site start',
                         'Site end', 'Site type', 'PCT', 'Conserved',
                         'Branch length score', 'UW BLS', 'Aligning species']

    site_info['UTR sequence'] = [x.replace('X', '')
                                 for x in UTRS_REF.loc[site_info['Gene ID']]
                                 ['UTR sequence']]
    site_info.loc[:, 'UTR BLS'] = list(BINS.loc[site_info['Gene ID']]['BLS'])

    # Truncate sites that overlap into the first 15 nts of the UTR
    new_sites = list(site_info['Site type'])
    new_starts = list(site_info['Site start'])
    new_pct = list(site_info['PCT'])
    new_conserved = list(site_info['Conserved'])
    for i in range(len(new_sites)):
        if new_starts[i] == options.RIBOSOME_SHADOW - 1:
            site = new_sites[i]
            if site == '8mer-1a':
                new_sites[i] = '7mer-1a'
                new_starts[i] += 1
            elif site == '7mer-m8':
                new_sites[i] = '6mer'
                new_starts[i] += 1
                new_pct[i] = 0
                new_conserved[i] = 0

    site_info['Site type'] = new_sites
    site_info['Site start'] = new_starts
    site_info['PCT'] = new_pct
    site_info['Conserved'] = new_conserved
    site_info = site_info[site_info['Site start'] >= options.RIBOSOME_SHADOW]

    print "{} seconds\n".format(time.time() - t0)

    # Write to file
    print "Writing outfile..."
    t0 = time.time()

    site_info.to_csv(options.OUT_FILE, sep='\t', index=False, float_format='%.6f')

    print "{} seconds\n".format(time.time() - t0)

    print "Total time for finding sites: {}\n".format(time.time() - T0)
