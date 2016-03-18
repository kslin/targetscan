import sys
import time

import concurrent.futures
import numpy as np
import pandas as pd

import config
import find_sites_helpers as fsh
import tasks


def main(SEED_FILE, BIN_FILE, UTR_FILE, OUT_FILE):
    T0 = time.time()

    # disable chained assignment warning
    pd.options.mode.chained_assignment = None

    # import seed information
    print "Importing seed, bin, and utr files..."
    t0 = time.time()

    SEED_TO_SPECIES, SEEDS = fsh.import_seeds(SEED_FILE)
    num_seeds = len(SEEDS)

    # Import gene to bin information
    BINS = pd.read_csv(BIN_FILE, sep='\t', header=None)
    BINS.columns = ['Gene', 'BLS', 'Bin']
    BINS = BINS.set_index('Gene')

    # import UTR information
    UTRS, UTRS_REF = fsh.import_utrs(UTR_FILE)

    print "{} seconds\n".format(time.time() - t0)


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
    site_df['Num sites'] = [fsh.occurrences(utr, fsh.rev_comp(seed[:-1]))
                            for (utr, seed)
                            in zip(site_df['UTR sequence'], site_df['Seed'])]
    site_df = site_df[site_df['Num sites'] > 0]

    # group the data by gene and pass to a helper function to get site info
    groups = site_df.groupby('Gene ID')
    num_genes = len(groups)
    site_info = []

    # if indicated by the user, use the parallel implementation
    if config.FUTURES:
        executor = concurrent.futures.ProcessPoolExecutor()
        futures = []
        for gene, group in groups:
            futures.append(executor.submit(tasks.get_all_site_info, gene,
                                           group, UTRS.loc[[gene]],
                                           SEED_TO_SPECIES))

            # add a sleep to prevent the executor from getting clogged
            time.sleep(0.0001)

        # as jobs are completed, add the information
        for future in concurrent.futures.as_completed(futures):
            site_info += future.result()

        executor.shutdown()

    # otherwise, use the nonparallel implementation
    else:
        for gene, group in groups:
            site_info += tasks.get_all_site_info(gene, group,
                                                 UTRS.loc[[gene]],
                                                 SEED_TO_SPECIES)

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
    print "{} seconds\n".format(time.time() - t0)


    # Write to file
    print "Writing outfile..."
    t0 = time.time()

    site_info.to_csv(OUT_FILE, sep='\t', index=False)

    print "{} seconds\n".format(time.time() - t0)

    print "Total time for finding sites: {}\n".format(time.time() - T0)


if __name__ == '__main__':

    # Get file of aligned UTRS and the file for writing output
    SEED_FILE, BIN_FILE, UTR_FILE, OUT_FILE = sys.argv[1:]

    # run main code
    main(SEED_FILE, BIN_FILE, UTR_FILE, OUT_FILE)
