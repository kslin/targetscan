import numpy as np
import pandas as pd


def calculate_TS7_scores(features_original, coeffs):
    
    # copy features dataframe and index by gene and miRNA
    features = features_original.copy()
    features.columns = [x.replace(' ','_') for x in features.columns]
    # features = features.set_index(keys=['Gene_ID', 'miRNA_family']).sort_index()

    # one-hot encode site-type
    ts7_stypes = list(features['Site_type'].unique())
    for stype in ts7_stypes:
        features[stype] = (features['Site_type'] == stype).astype(float)

    # set upper bound for scores
    upper_bound_dict = {
        '8mer-1a': -0.03,
        '7mer-m8': -0.02,
        '7mer-1a': -0.01,
        '6mer': 0.0
    }

    features['upper_bound'] = [upper_bound_dict[x] for x in features['Site_type']]

    # normalize features and weight with coefficients
    weighted_features = []
    for stype, group in features.groupby('Site_type'):
        group['score'] = 0
        for row in coeffs.iterrows():
            param = row[0]
            if param == 'Intercept':
                group['Intercept'] = row[1]['{} coeff'.format(stype)]
            else:
                feat_min, feat_max, feat_coeff = row[1][['{} min'.format(stype), '{} max'.format(stype), '{} coeff'.format(stype)]]
                vals = group[param]
                group[param] = feat_coeff * (vals - feat_min) / (feat_max - feat_min)
        weighted_features.append(group)
    weighted_features = pd.concat(weighted_features)
    weighted_features = weighted_features[['Gene_ID', 'miRNA_family', 'upper_bound', 'AIR'] + list(coeffs.index)]

    # add up all weighted features to get a raw score
    weighted_features['score'] = np.sum(weighted_features[list(coeffs.index)].values, axis=1)
    weighted_features['weighted_score'] = weighted_features['AIR'] * weighted_features['score']

    # bound score according to Agarwal 2015
    weighted_features['bounded_score'] = np.minimum(weighted_features['upper_bound'], weighted_features['score'])
    weighted_features['weighted_bounded_score'] = weighted_features['AIR'] * weighted_features['bounded_score']

    # aggregate scores for each transcript, miRNA pair
    pred_df = weighted_features.groupby(['Gene_ID', 'miRNA_family']).agg({
        'score': np.sum,
        'weighted_score': np.sum,
        'bounded_score': np.sum,
        'weighted_bounded_score': np.sum
    })

    return weighted_features, pred_df
