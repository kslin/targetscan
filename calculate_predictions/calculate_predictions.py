#!/usr/bin/python2

from optparse import OptionParser
import os
import time

import numpy as np
import pandas as pd

from tasks import calculate_TS7_scores


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--feature_file", dest="FEATURE_FILE", help="output from calculate_features")
    parser.add_option("--coeff_file", dest="COEFF_FILE", help="file with model coefficients")
    parser.add_option("--out_dir", dest="OUT_DIR", help="prefix for writing outputs")

    (options, args) = parser.parse_args()

    T0 = time.time()

    # disable chained assignment warning
    pd.options.mode.chained_assignment = None

    # multiply by parameters published in Agarwal 2015 and write final predictions
    print "Calculating context++ scores..."

    FEATURES = pd.read_csv(options.FEATURE_FILE, sep='\t')
    COEFFS = pd.read_csv(options.COEFF_FILE, sep='\t', index_col='Feature')

    WEIGHTED_FEATURES, PRED_DF = calculate_TS7_scores(FEATURES, COEFFS)

    WEIGHTED_FEATURES.to_csv(os.path.join(options.OUT_DIR, 'weighted_features.txt'), sep='\t', index=False, float_format='%.6f')
    PRED_DF.to_csv(os.path.join(options.OUT_DIR, 'predictions.txt'), sep='\t', float_format='%.6f')

    print '{} seconds\n'.format(time.time() - T0)
