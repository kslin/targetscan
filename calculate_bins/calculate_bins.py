import csv
import sys
import time

from Bio import Phylo
import concurrent.futures
import pandas as pd
# from tqdm import tqdm

import config
from tasks import get_median_bls


T0 = time.time()

UTR_FILE = sys.argv[1]
OUT_FILE = sys.argv[2]

sys.stdout.write('Processing UTRS... ')

UTRS = pd.read_csv(UTR_FILE,sep='\t',header=None)
UTRS.columns = ['Gene ID','Species ID','UTR sequence']
UTRS[['Gene ID','Species ID','UTR sequence']] = UTRS[['Gene ID','Species ID','UTR sequence']].astype(str)
UTRS = UTRS.groupby('Gene ID')

sys.stdout.write('{} seconds\n'.format(time.time()-T0))
T0 = time.time()
sys.stdout.write('Assigning bins... ')
bins = []
num_utrs = len(UTRS)
if config.FUTURES:
    executor = concurrent.futures.ProcessPoolExecutor()
    # futures = [executor.submit(get_median_bls, i, gene, utr_group) for i,(gene, utr_group) in enumerate(UTRS)]
    futures = []
    for i,(gene,utr_group) in enumerate(UTRS):
        futures.append(executor.submit(get_median_bls, i, gene, utr_group))
        if (i % 1000) == 0:
            print '{}/{}'.format(i,num_utrs)

    # for future in tqdm(concurrent.futures.as_completed(futures), total=len(UTRS), unit='UTR', leave=True):
    for future in concurrent.futures.as_completed(futures):
        bins.append(future.result())

    executor.shutdown()
else:
    # for gene, utr_group in tqdm(UTRS, unit='UTR', leave=True):
    for i,(gene, utr_group) in enumerate(UTRS):
        bins.append(get_median_bls(i, gene, utr_group))

# sys.stdout.write(bins
sys.stdout.write('{} seconds\n'.format(time.time()-T0))
# mydf = pd.DataFrame(bins)

sys.stdout.write('Writing to file... ')
T0 = time.time()
# sys.stdout.write(bins
# mydf.to_csv(OUT_FILE,sep='\t',index=False,header=None)
with open(OUT_FILE,'wb') as f:
    writer = csv.writer(f,delimiter='\t')
    writer.writerows(bins)
sys.stdout.write('{} seconds\n'.format(time.time()-T0))


