import random
import sys

import pandas as pd

import linear_regression_helpers

training = pd.read_csv('../../data/training.txt',sep='\t')
testing = pd.read_csv('../../data/testing.txt',sep='\t')

training_one_site = training[training['Num sites'] == 1]

site_specific = linear_regression_helpers.SiteSpecific()
site_specific.fit(training_one_site, cols)

mtraining = linear_regression_helpers.organize0(training,site_specific)
multisite = linear_regression_helpers.MultiSite()
multisite.fit(mtraining)
mtraining['multiscore'] = multisite.predict(mtraining)
print multisite.get_coefs()

straining = linear_regression_helpers.prepare_data_for_site_specific(mtraining,training)
for i in range(10):
    site_specific = SiteSpecific()
    site_specific.fit(straining,cols)
    
    mtraining = linear_regression_helpers.organize0(training,site_specific)
    multisite = linear_regression_helpers.MultiSite()
    multisite.fit(mtraining)
    mtraining['multiscore'] = multisite.predict(mtraining)
    print multisite.get_coefs()

    straining = linear_regression_helpers.prepare_data_for_site_specific(mtraining,training)



mtesting = linear_regression_helpers.organize0(testing,site_specific)
print np.nanmean(mtesting['logFC'])
predicted = np.add(multisite.predict(mtesting),mtesting['Single site score'].values)
print get_r_squared(list(mtesting['Single site score']),list(mtesting['logFC']))
print get_r_squared(predicted,mtesting['logFC'])
mtesting['combined score'] = predicted

mtesting.to_csv('../../data/model_results.txt',sep='\t',index=False)

