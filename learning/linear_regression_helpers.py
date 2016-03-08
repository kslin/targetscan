import copy
import random
import time
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression,Lasso
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor



FEATURE_FILE = '../../data/big_output.txt'
LOGFC_FILE1 = '../../data/Supplementary1.csv'
LOGFC_FILE2 = '../../data/Supplementary2.csv'
GENE_FILE = '../../targetscan_files/Gene_info.txt'
SEED_FILE = '../../seed_dict.csv'
TS7_human_file = '../../Summary_Counts_Human7.txt'

cols = ['Threep score','Local AU score','Min dist score','TA', 'SPS',
        'UTR length score','Off6m score','SA',
        'siRNA 1A', 'siRNA 1C', 'siRNA 1G', 'siRNA 8A',
        'siRNA 8C', 'siRNA 8G', 'site 8A', 'site 8C', 'site 8G', 
        'PCT','ORF length','ORF 8mers']


def get_r_squared(x,y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return r_value**2

def get_logfc(df,gene,seed):
    if gene not in df.index:
        return float('NaN')
    if seed not in df.columns:
        return float('NaN')
    return float(df.loc[gene][seed])

def get_num_coop(locs,types):
    upper = 46
    lower = 13
    locs = [l for (l,t) in zip(locs,types) if t != '6mer']
    if len(locs) < 2:
        return 0
    num_coop = 0
    for i in range(len(locs)-1):
        diff = locs[i+1] - locs[i]
        if (diff < 46) & (diff > 13):
            num_coop += 1
    return num_coop

def organize0(data,site_specific):
        subdf = copy.deepcopy(data)
        subdf['single site score'] = site_specific.predict(data,cols)
        subdf = subdf.groupby('ID').agg({'single site score': lambda x: tuple(x),
                                                       'logFC': lambda x: tuple(x)[0],
                                                       'Site type': lambda x: tuple(x),
                                                       'Site end': lambda x: tuple(x),
                                                       'Num sites': lambda x: tuple(x)[0],
                                                       'UTR length score': np.nanmean,
                                                       'special id':lambda x: tuple(x)})
        subdf['Num coop'] = [get_num_coop(locs,sites) for (locs,sites) in zip(subdf['Site end'],subdf['Site type'])]
        subdf['Num sites score'] = np.divide(subdf['Num sites'],subdf['UTR length score'])
        subdf['Num coop score'] = np.divide(subdf['Num coop'],subdf['UTR length score'])
        subdf['Single site score'] = [sum(x) for x in subdf['single site score']]
        return subdf

def organize(data,site_specific,old_subdf):
        subdf = copy.deepcopy(data)
        subdf['single site score'] = site_specific.predict(data,cols)
        subdf = subdf[['ID','single site score']]
        subdf = subdf.groupby('ID').agg(lambda x: tuple(x))

        old_len = len(subdf)
        subdf = pd.concat([subdf,old_subdf.drop('single site score',1)],axis=1,join='inner')

        assert len(subdf) == old_len, '{},{}'.format(len(subdf),old_len)

        subdf['Single site score'] = [sum(x) for x in subdf['single site score']]

        return subdf

def prepare_data_for_site_specific(data, original):
    t0 = time.time()
    adjusted_logFC = list(np.subtract(data['logFC'],data['multiscore']))
    zipped = zip(data['Single site score'],data['single site score'],adjusted_logFC)
    new_logfcs = [[(x/site_score)*logfc if (site_score != 0) else (logfc/len(site_scores)) for x in site_scores] for (site_score,site_scores,logfc) in zipped]
    # t0 = time.time()
    ids = list(sum(list(data['special id']), ()))
    # print 'woah4 {}'.format(time.time()-t0)
    # t0 = time.time()
    new_logfcs = np.array(sum(new_logfcs,[]))
    # print 'woah5 {}'.format(time.time()-t0)
    # t0 = time.time()

    return_df = copy.deepcopy(original)
    return_df['logFC'] = new_logfcs[np.argsort(ids)]

    return return_df



#### Training methods #####

class SiteSpecific(object):

    def __init__(self, normalize=False):
        self.ols_models = [LinearRegression(normalize=normalize) for i in range(4)]
        self.site_maxes = [-0.03,-0.02,-0.01,0]

    def fit(self, data, cols):
        for i,sitetype in enumerate(['8mer-1a', '7mer-m8', '7mer-1a', '6mer']):
            subdf = data[data['Site type'] == sitetype]
            if len(subdf) == 0:
                continue
            self.ols_models[i].fit(subdf[cols],subdf['logFC'].values)

    def predict(self, data, cols):
        mydata = copy.deepcopy(data)
        mydata['ix'] = range(len(mydata))
        ixs, values = [], []
        for i,sitetype in enumerate(['8mer-1a', '7mer-m8', '7mer-1a', '6mer']):
            subdf = mydata[mydata['Site type'] == sitetype]
            if len(subdf) == 0:
                continue
            ixs += list(subdf['ix'])
            predictions = self.ols_models[i].predict(subdf[cols])
            values += [min(self.site_maxes[i],p) for p in predictions]

        return np.array(values)[np.argsort(ixs)]

    def score(self, data, cols):
        actual = list(data['logFC'])
        predicted = list(self.predict(data, cols))
        scores = []
        for i,sitetype in enumerate(['8mer-1a', '7mer-m8', '7mer-1a', '6mer']):
            subdf = data[data['Site type'] == sitetype]
            if len(subdf) == 0:
                continue
            scores.append(self.ols_models[i].score(subdf[cols].values,subdf['logFC'].values))
        
        scores.append(get_r_squared(actual, predicted))

        return scores
    
    def get_coefs(self):
        coefs = []
        intercepts = []
        for i in range(4):
            try:
                coefs.append(list(self.ols_models[i].coef_))
                intercepts.append(self.ols_models[i].intercept_)
            except:
                continue
        
        return coefs,intercepts


class MultiSite(object):
    
    def __init__(self, normalize=False):
        self.ols_model = LinearRegression(normalize=normalize)
    
    def fit(self,data):
        Ytrain = np.subtract(data['logFC'].values,data['Single site score'].values)
        self.ols_model.fit(data[['Num sites score','Num coop score']],Ytrain)
        predicted = self.ols_model.predict(data[['Num sites score','Num coop score']])
        print 'Training score: {}'.format(get_r_squared(np.add(predicted,data['Single site score']),data['logFC']))
    
    def predict(self,data):
        Ytest = np.subtract(data['logFC'].values,data['Single site score'].values)
        return self.ols_model.predict(data[['Num sites score','Num coop score']])
    
    def score(self,data):
        predicted = np.add(self.predict(data),data['Single site score'].values)
        print get_r_squared(list(data['Single site score']),list(data['logFC']))
        return get_r_squared(predicted,data['logFC'])

    def get_coefs(self):
        return self.ols_model.coef_, self.ols_model.intercept_







