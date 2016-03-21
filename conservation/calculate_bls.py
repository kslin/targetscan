import csv
import string
import sys
import time

from Bio import Phylo
import concurrent.futures
import pandas as pd

import config
from tasks import get_bls


def main(UTR_FILE, OUT_FILE):
    T0 = time.time()

    # PROCESS UTR FILE
    print 'Processing UTRS... '
    t0 = time.time()

    # Read UTR file into a pandas dataframe and group UTRs by gene
    UTRS = pd.read_csv(UTR_FILE, sep='\t', header=None)
    UTRS.columns = ['Gene ID', 'Species ID', 'UTR sequence']
    UTRS[['Gene ID', 'Species ID', 'UTR sequence']] \
        = UTRS[['Gene ID', 'Species ID', 'UTR sequence']].astype(str)
    UTRS = UTRS.groupby('Gene ID')

    print '{} seconds\n'.format(time.time()-t0)

    # Calculate branch length scores
    print 'Calculating bls... '
    t0 = time.time()
    bins = []
    num_utrs = len(UTRS)

    # Use the parallel version of the code if indicated by the user
    if config.FUTURES:
        executor = concurrent.futures.ProcessPoolExecutor()
        futures = []

        # Spin off jobs to be done in parallel
        for i, (gene, utr_group) in enumerate(UTRS):
            futures.append(executor.submit(get_bls, gene, utr_group))

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
            bins.append(get_bls(gene, utr_group))

    print '{} seconds\n'.format(time.time()-t0)

    # WRITE RESULTS TO A FILE
    print 'Writing to file... '
    t0 = time.time()

    with open(OUT_FILE, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(bins)
    print '{} seconds\n'.format(time.time()-t0)

    print 'Total time for calculating bls: {}\n'.format(time.time()-T0)


if __name__ == '__main__':

    # Get file of aligned UTRS and the file for writing output
    UTR_FILE = sys.argv[1]
    OUT_FILE = sys.argv[2]

    # run main code
    main(UTR_FILE, OUT_FILE)
