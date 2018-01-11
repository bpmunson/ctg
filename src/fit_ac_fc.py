"""Initializes abundance and fitnesses for the iterative least squares

This module initializes abundance and fitness for iterative least squares later.

TODO:
    * Update docstrings
    * Separate conversion tools (to change existing data structure to new ones)
    * Change variable names to something more descriptive

"""

import os
import numpy as np
from numpy.matlib import repmat
import pandas as pd
from scipy.stats import t

def fit_ac_fc(ab, freq):
    pass

def load_abundance_thresholds(fn, sep='\s+'):
    """Loads abundance thresholds.

    """

    ab = pd.read_csv(fn, sep=sep, index_col=0)

    return ab

def load_timepoint_counts(fn, sep="\s+"):
    """Loads timepoint counts.

    """

    tps = pd.read_csv(fn, sep=sep)

    return tps

def convert_abundance_thresholds(df, sep="_"):
    df.index = pd.MultiIndex.from_tuples(
        [(int(i.split('_')[2]), int(i.split('_')[1][1:])) for i in df.index],
        names=['reps', 'time'])

    return df.T

def convert_timepoint_counts(df):
    data = df.iloc[:, 5:]
    names = df.iloc[:,:5]

    idx = pd.MultiIndex.from_tuples(
        [(int(i.split('_')[2]), int(i.split('_')[1][1:])) for i in data.columns],
        names=['reps', 'time'])

    data.columns = idx

    return data, names

def _cov(x,y, axis=0):
    return np.ma.mean(x*y, axis=axis) - (np.ma.mean(x, axis=axis)*np.ma.mean(y, axis=axis))

def fit_ac_fc_np(counts, ab, times, n_good=2):
    #TODO: Everything
    #TODO: Missing local_fdr

    bad = (counts > ab).sum(axis=2) < n_good
    allbad = bad.all(axis=0)

    counts_masked = np.ma.array(data=counts, mask=~(counts > ab))
    time_masked = np.ma.array(data=times, mask=~(counts > ab))

    mean_counts = np.ma.mean(counts_masked, axis=2).data
    mean_time = np.ma.mean(time_masked, axis=2).data
    var_time = np.ma.var(time_masked,axis=2).data
    f = _cov(counts_masked, time_masked, axis=2).data

    fc = np.divide(f.sum(axis=0), var_time.sum(axis=0))
    fc[allbad] = 0

    ac = mean_counts - (fc*mean_time)
    ac[bad] = counts[bad,0]

    alpha = -np.log2(np.power(2,ac).sum(axis=1))[...,np.newaxis]
    ac = ac + alpha

    lmbda = -np.log2(np.power(2, ac[...,np.newaxis] + fc[np.newaxis,...,np.newaxis]*times).sum(axis=1))
    xfit = ac[...,np.newaxis] + fc[np.newaxis,...,np.newaxis]*times + lmbda[...,np.newaxis].transpose(0,2,1)

    df = (counts > ab).sum(axis=2).sum(axis=0) - 2
    numerator = np.sqrt(np.power(xfit - counts_masked, 2).sum(axis=2).sum(axis=0))
    denominator = np.sqrt(np.power(time_masked - time_masked.mean(axis=2)[...,np.newaxis], 2).sum(axis=2).sum(axis=0))

    sdfc = np.divide(numerator,denominator)

    df_masked = np.ma.array(data=df, mask=~(df > 0))

    tstat_masked = np.ma.divide(fc, np.ma.divide(sdfc, np.ma.sqrt(df_masked)))
    tstat_masked.data[tstat_masked.mask] = 1

    tstat = tstat_masked.data
    p_t = np.array([2*t.cdf(i,j) if j > 0 else 1 for i,j in zip(tstat, df)])

    return ac, fc, sdfc, p_t, df, allbad

if __name__ == "__main__":
    from config import config

    ab = load_abundance_thresholds(os.path.join(config.A549_test, \
            "A549_abundance_thresholds.txt"))
    tps = load_timepoint_counts(os.path.join(config.A549_test, \
            "A549_timepoint_counts.txt"))

    _tps, names = convert_timepoint_counts(tps)
    _abundance = convert_abundance_thresholds(ab)

    #For Testing only
    _tps = _tps.iloc[:100,:]
    _abundance = _abundance.iloc[:100,:]


    times = np.array([[3, 14, 21, 28], [3,14,21,28]])
    #times = np.transpose(times[...,np.newaxis], axes=[0,2,1])
    times = repmat(times, 100, 1).reshape((2,100,4))

    fit_ac_fc_np(np.array([_tps[1].as_matrix(), _tps[2].as_matrix()]),
                 np.array([_abundance[1].as_matrix(), _abundance[2].as_matrix()]),
                 times)
