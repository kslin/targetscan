import copy
import random

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


def get_r_squared(x,y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return r_value**2

def get_logfc(df,gene,seed):
    if gene not in df.index:
        return float('NaN')
    if seed not in df.columns:
        return float('NaN')
    return float(df.loc[gene][seed])

def get_data():
    GENE_INFO = pd.read_csv(GENE_FILE,sep='\t').drop(['Gene description','Species ID'],1)
    GENE_INFO = GENE_INFO.groupby('Gene symbol').agg(lambda x:tuple(x))
    GENE_INFO.loc[:,'Isoform ratio'] = [float(max(x))/np.nansum(x) for x in GENE_INFO['3P-seq tags + 5']]
    GENE_INFO.loc[:,'Transcript ID'] = [x[y.index(1)] for (x,y) in zip(GENE_INFO['Transcript ID'],GENE_INFO['Representative transcript?'])]
    GENE_INFO = GENE_INFO[['Transcript ID','Isoform ratio']]

    SEED_INFO = pd.read_csv(SEED_FILE,sep='\t')
    SEED_DICT = {}
    for row in SEED_INFO.iterrows():
        SEED_DICT[row[1]['col']] = row[1]['seed']

    logFCs1 = pd.read_csv(LOGFC_FILE1)
    logFCs1['Transcript ID'] = list(GENE_INFO.loc[logFCs1['Gene symbol']]['Transcript ID'])
    logFCs1['Isoform ratio1'] = list(GENE_INFO.loc[logFCs1['Gene symbol']]['Isoform ratio'])
    logFCs1 = logFCs1.drop(['Used in training','RefSeq ID','Gene symbol'],1).dropna(subset=['Transcript ID'])
    logFCs1 = logFCs1.drop_duplicates(subset=['Transcript ID']).set_index('Transcript ID')
    logFCs1.columns = [SEED_DICT[x] if x in SEED_DICT else x for x in logFCs1.columns]

    logFCs2 = pd.read_csv(LOGFC_FILE2)
    logFCs2['Transcript ID'] = list(GENE_INFO.loc[logFCs2['Gene symbol']]['Transcript ID'])
    logFCs2['Isoform ratio2'] = list(GENE_INFO.loc[logFCs2['Gene symbol']]['Isoform ratio'])
    logFCs2 = logFCs2.drop(['RefSeq ID','Gene symbol'],1).dropna(subset=['Transcript ID'])
    logFCs2 = logFCs2.drop_duplicates(subset=['Transcript ID']).set_index('Transcript ID')
    logFCs2.columns = [SEED_DICT[x] if x in SEED_DICT else x for x in logFCs2.columns]

    logFCs = pd.concat([logFCs1,logFCs2],axis=1,join='outer')


    features = pd.read_csv(FEATURE_FILE, sep='\t')

    features['logFC'] = [get_logfc(logFCs,gene,seed) for (gene,seed) in zip(features['Gene ID'],features['Seed'])]
    features = features.dropna()

    isoform_ratios = [y if np.isnan(x) else x for (x,y) in zip(logFCs['Isoform ratio1'],logFCs['Isoform ratio2'])]
    logFCs['Isoform ratio'] = isoform_ratios

    features['Isoform ratio'] = list(logFCs.loc[features['Gene ID']]['Isoform ratio'])
    features['Num sites'] = [1]*len(features)
    num_sites = features.groupby(['Gene ID','Seed']).agg({'Num sites':np.nansum})
    features['Num sites'] = list(num_sites.loc[zip(features['Gene ID'],features['Seed'])]['Num sites'])

    return features

def get_ts7():
    targetscan7 = pd.read_csv(TS7_human_file,sep='\t',usecols=['Transcript ID','Gene Symbol','miRNA family',
                                                                                                    'Total num conserved sites',
                                                                                                    'Cumulative weighted context++ score'])
    targetscan7 = targetscan7[targetscan7['Transcript ID'].isin(GENE_INFO['Transcript ID'])]
    targetscan7 = targetscan7[targetscan7['miRNA family'].isin(all_seeds)]
    targetscan7 = targetscan7.rename(columns={'Gene Symbol':'Gene symbol'})
    targetscan7.loc[:,'logFC'] = [get_logfc())(targetscan7['Gene symbol'],targetscan7['miRNA family'])
    targetscan7 = targetscan7.dropna(subset=['logFC'])
    targetscan7 = targetscan7.set_index('miRNA family')
    
    return targetscan7

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





def organize0(data,site_specific):
        subdf = copy.deepcopy(data)
        subdf['special id'] = range(len(data))
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

def organize(data,site_specific):
        subdf = copy.deepcopy(data)
        subdf['special id'] = range(len(data))
        subdf['single site score'] = site_specific.predict(data,cols)
        subdf = subdf.groupby(['Gene ID','Seed']).agg({'single site score': lambda x: tuple(x),
                                                       'logFC': np.nansum,
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

def prepare_data_for_site_specific(data, original):
    data = copy.deepcopy(data)
    data['adjusted logFC'] = np.subtract(data['logFC'],data['multiscore'])
    new_logfcs = [[(x/site_score)*logfc if (site_score != 0) else (logfc/len(site_scores)) for x in site_scores] for (site_score,site_scores,logfc) in zip(data['Single site score'],
                                                                                                                                 data['single site score'],
                                                                                                                                 data['adjusted logFC'])]
    ids = reduce(lambda x,y: x+y, list(data['special id']))
    new_logfcs = np.array(reduce(lambda x,y: x+y, new_logfcs))
    return_df = copy.deepcopy(original)
    return_df['logFC'] = new_logfcs[np.argsort(ids)]

    return return_df






