import copy
import random
import sys
import time

import numpy as np
import pandas as pd

import linear_regression_helpers

TRAIN_FILE = '../../learning_files/training.txt'
TEST_FILE = '../../learning_files/testing.txt'
OUT_FILE = '../../learning_files/model_results.txt'


cols = ['Threep score','Local AU score','Min dist score','TA', 'SPS',
        'UTR length score','Off6m score','SA',
        'siRNA 1A', 'siRNA 1C', 'siRNA 1G', 'siRNA 8A',
        'siRNA 8C', 'siRNA 8G', 'site 8A', 'site 8C', 'site 8G', 
        'PCT','ORF length','ORF 8mers']


training = pd.read_csv(TRAIN_FILE,sep='\t')
training['special id'] = range(len(training))
training = training.iloc[0:10000]

testing = pd.read_csv(TEST_FILE,sep='\t')
testing['special id'] = range(len(testing))

training_one_site = training[training['Num sites'] == 1]

originalt = training.groupby(['Gene ID','Seed']).agg({'single site score': np.nansum,'logFC':lambda x: tuple(x)[0]})
print np.nanmean(originalt['logFC'])

site_specific = linear_regression_helpers.SiteSpecific()
site_specific.fit(training_one_site, cols)

mtraining0 = linear_regression_helpers.organize0(training,site_specific)
mtraining = copy.deepcopy(mtraining0)
multisite = linear_regression_helpers.MultiSite()
multisite.fit(mtraining)
mtraining['multiscore'] = multisite.predict(mtraining)
print multisite.get_coefs()

straining = linear_regression_helpers.prepare_data_for_site_specific(mtraining,training)
t0 = time.time()
for i in range(10):
    site_specific = linear_regression_helpers.SiteSpecific()
    site_specific.fit(straining,cols)

    mtraining = linear_regression_helpers.organize(training,site_specific,mtraining0)

    multisite = linear_regression_helpers.MultiSite()

    multisite.fit(mtraining)
    mtraining['multiscore'] = multisite.predict(mtraining)

    print multisite.get_coefs()

    straining = linear_regression_helpers.prepare_data_for_site_specific(mtraining,training)




mtesting = linear_regression_helpers.organize0(testing,site_specific)
print np.nanmean(mtesting['logFC'])
predicted = np.add(multisite.predict(mtesting),mtesting['Single site score'].values)
print linear_regression_helpers.get_r_squared(list(mtesting['Single site score']),list(mtesting['logFC']))
print linear_regression_helpers.get_r_squared(predicted,mtesting['logFC'])
mtesting['combined score'] = predicted

# mtesting.to_csv(OUT_FILE,sep='\t')

