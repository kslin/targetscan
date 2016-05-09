# import matplotlib.pyplot as plt
import collections
import numpy as np
import pandas as pd
import pickle
import pystan
import scipy
from sklearn.linear_model import LinearRegression
import sys

infile = '../../learning_files/features_norm_1site.txt'
feat_names = ['Threep score', 'Local AU score', 'Min dist score', 'UTR length score',
              'Off6m score', 'SA',
              # 'siRNA 1A', 'siRNA 1C', 'siRNA 1G',
              # 'siRNA 8A', 'siRNA 8C', 'siRNA 8G',
              # 'site 8A', 'site 8C', 'site 8G', 
              'PCT', 
              # 'Branch length score', 'UTR BLS',
              'TA', 'SPS', 'ORF length', 'ORF 8mers']
# feat_names = ['Threep score', 'Local AU score', 'Min dist score', 'UTR length score']

features = pd.read_csv(infile, sep='\t')
features = features[features['Site type'] == '8mer-1a']

np.random.seed(1)
all_indices = set(range(len(features)))
testing_indices = list(np.random.choice(range(len(features)), 5000, replace=False))
training_indices = list(all_indices - set(testing_indices))
training = features.iloc[training_indices]
testing = features.iloc[testing_indices]
# features = features.iloc[np.random.choice(range(len(features)), 2000)]
train_experiments = list(set(list(training['Seed'])))
train_exp_dict = {seed:(i+1) for i,seed in enumerate(train_experiments)}
training.loc[:,'group'] = [train_exp_dict[seed] for seed in training['Seed']]

test_experiments = list(set(list(testing['Seed'])))
test_exp_dict = {seed:(i+1) for i,seed in enumerate(test_experiments)}
testing.loc[:,'group'] = [test_exp_dict[seed] for seed in testing['Seed']]

Xtrain = training[feat_names]
Xtrain.loc[:,'intercept'] = [1]*len(Xtrain)
Xtrain = Xtrain.values
Ytrain = training['-logFC'].values
train_groups = training['group'].values

Xtest = testing[feat_names]
Xtest.loc[:,'intercept'] = [1]*len(Xtest)
Xtest = Xtest.values
Ytest = testing['-logFC'].values
test_groups = testing['group'].values

# Xtest = np.random.random((100, 5))
# true_coefs = np.array([1,2,3,4,5]).T
# noise = np.random.random(100)
# Ytest = np.dot(Xtest, true_coefs) + noise
# print Xtest.shape
# print Ytest.shape
# test_experiments = [1]
# test_groups = [1]*100

# coefs1 = np.array([0.06, 0.16, -0.11, 0.3, -0.13, 0.2, -0.03, -0.06, -0.02, 6.2e-3,
#                    0.1, 0.05, 1.38, -0.14, 0.1, 0.2, 0.21, 0.21, -0.1, -0.69])

# coefs2 = np.array([0.32, 1.45, 0.19, 1.36, -0.62, 1.46, -0.02, 2.7e-3, -0.47, -0.11,
#                    0.3, 0.2, 1.2, 1.62, 1.13, 1.03, 3.01, 1.75, -0.49, -7.37])

# coefs1 = np.array([0.56, 0.58, -0.01, 0.21, -0.23, -0.05, 0.09, -0.07, -0.07, -0.2,
#                    0.06, 0.05, -2.65, 0.19, 0.03, 0.13, -0.07, -0.03, -0.09, 0.03])

# coefs2 = np.array([0.55, 0.95, 0.22, 0.25, -0.2, 0.01, -7.0e-4, -0.03, -0.17, 0.09,
#                    0.32, 0.07, -11.23, 1.16, 0.7, 0.35, 0.19, 0.69, -0.28, -1.61 ])

# predicted_y = np.multiply(np.dot(Xtest, coefs1), np.dot(Xtest, coefs2))
# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Ytest, predicted_y)
# print r_value**2
# sys.exit()

# print np.nansum(X, axis=0)
# sys.exit()

model = LinearRegression()
model.fit(Xtrain, Ytrain)
print 'OLS score: {}'.format(model.score(Xtest, Ytest))

# traces = np.load('pystan_trace_effs.npy')
# A = str(traces)
# B = A.split(')')[0]
# print B.split(',')
# coefs = traces['coefs']
# predicted_y = np.nanmean(np.dot(coefs, Xtest.T), axis=0)
# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Ytest, predicted_y)
# print r_value**2
sys.exit()

# print X.shape
# print Y.shape
# sys.exit()

# print X
# print np.sum(X, axis=0)
# sys.exit()

code = """
data {
    int<lower=0> Ndata;
    int<lower=0> Nfeats;
    int<lower=0> Ngroups;
    matrix[Ndata, Nfeats] X;
    vector[Ndata] Y;
    int<lower=1,upper=Ngroups> group_ids[Ndata];
}
parameters {
    vector[Nfeats] coefs;
    real<lower=0> stdev;
    vector<lower=0,upper=1>[Ngroups] effs;
}
transformed parameters {
    vector[Ndata] mu;
    vector<lower=0, upper=1>[Ndata] eff_list;
    for (n in 1:Ndata){
        eff_list[n] <- effs[group_ids[n]];
    }
    mu <- (X*coefs) .* eff_list;
}
model {
    effs ~ beta(1,1);
    Y ~ normal(mu, stdev);
}
"""

data = {'Ndata': len(Xtrain),
        'Nfeats': len(feat_names)+1,
        'Ngroups': len(train_experiments),
        'X': Xtrain,
        'Y': Ytrain,
        'group_ids': train_groups}

fit = pystan.stan(model_code=code, data=data)
# print fit
traces = fit.extract()
np.save("pystan_trace_coefs8.npy", traces['coefs'])
np.save("pystan_trace_effs8.npy", traces['effs'])
np.save("pystan_trace_stdev8.npy", traces['stdev'])
coefs = traces['coefs']

predicted_y = np.nanmean(np.dot(coefs, Xtest.T), axis=0)
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Ytest, predicted_y)
print r_value**2

print fit

# predicted_y = np.dot(Xtest, coefs)
# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Ytest, predicted_y)
# print r_value**2

# sm1 = pystan.StanModel(model_code=code1)
# fit1 = sm1.sampling(data=data)
# print(fit1)

# sm = pystan.StanModel(model_code=code)
# fit = sm.sampling(data=data)
# print(fit)

# # save it to the file 'model.pkl' for later use
# with open('model_short.pkl', 'wb') as f:
#     pickle.dump(sm, f)
