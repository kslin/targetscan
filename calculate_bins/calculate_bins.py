#!/usr/bin/python2

import csv
from optparse import OptionParser
import string
import sys
import time

from Bio import Phylo
import concurrent.futures
import pandas as pd

import bin_helpers
from tasks import get_median_bls


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--utr3_file", dest="UTR3_FILE", help="3\'UTR sequences in tab delimited format")
    parser.add_option("--tree_file", dest="TREE_FILE", help="phylogenetic tree in newick format")
    parser.add_option("--out", dest="OUT_FILE", help="where to write bin output")
    parser.add_option("--ref_species", dest="REF_SPECIES", help="reference species", default='9606')
    parser.add_option("--futures", dest="FUTURES", help="if true, run in parallel", default=False, action='store_true')

    (options, args) = parser.parse_args()

    T0 = time.time()

    # Process phylogenetic tree information
    SPECIES_TO_PATH = bin_helpers.parse_tree(options.TREE_FILE, options.REF_SPECIES)

    # Process UTR FILE
    print('Processing UTRS... ')
    t0 = time.time()

    # Read UTR file into a pandas dataframe and group UTRs by gene
    UTRS = pd.read_csv(options.UTR3_FILE, sep='\t', header=None)
    UTRS.columns = ['Gene ID', 'Species ID', 'UTR sequence']
    UTRS[['Gene ID', 'Species ID', 'UTR sequence']] \
        = UTRS[['Gene ID', 'Species ID', 'UTR sequence']].astype(str)
    UTRS = UTRS.groupby('Gene ID')

    print('{} seconds\n'.format(time.time()-t0))

    # ASSIGN UTRS TO BINS
    print('Assigning bins... ')
    t0 = time.time()
    bins = []
    num_utrs = len(UTRS)

    # Use the parallel version of the code if indicated by the user
    if options.FUTURES:
        executor = concurrent.futures.ProcessPoolExecutor()
        futures = []

        # Spin off jobs to be done in parallel
        for i, (gene, utr_group) in enumerate(UTRS):
            futures.append(executor.submit(get_median_bls, gene, utr_group, SPECIES_TO_PATH, options.REF_SPECIES))

            # Add sleep to prevent the executor from being overwhelmed
            time.sleep(0.0001)

        # Collect job results as they come in
        # since we do not care about the order in which jobs finish
        for future in concurrent.futures.as_completed(futures):
            bins.append(future.result())

        executor.shutdown()

    # Otherwise, use the non-parallel version
    else:
        for i, (gene, utr_group) in enumerate(UTRS):
            bins.append(get_median_bls(gene, utr_group, SPECIES_TO_PATH, options.REF_SPECIES))

    print('{} seconds\n'.format(time.time()-t0))

    # WRITE RESULTS TO A FILE
    print('Writing to file... ')
    t0 = time.time()

    with open(options.OUT_FILE, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(bins)
    print('{} seconds\n'.format(time.time()-t0))

    print('Total time for calculating bins: {}\n'.format(time.time()-T0))
