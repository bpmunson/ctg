import os
import sys

import numpy as np
from numpy.matlib import repmat
import pandas as pd
from scipy.stats import t

from config import config

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri

def _load_abundance_thresholds(fn, sep='\s+'):
    """Loads abundance thresholds.

    """

    ab = pd.read_csv(fn, sep=sep, index_col=0)

    return ab

def _load_timepoint_counts(fn, sep="\s+"):
    """Loads timepoint counts.

    """

    tps = pd.read_csv(fn, sep=sep)

    return tps

def _convert_abundance_thresholds(df, sep="_"):
    df.index = pd.MultiIndex.from_tuples(
        [(int(i.split('_')[2]), int(i.split('_')[1][1:])) for i in df.index],
        names=['reps', 'time'])

    return df.T

def _convert_timepoint_counts(df):
    data = df.iloc[:, 5:]
    names = df.iloc[:,:5]

    idx = pd.MultiIndex.from_tuples(
        [(int(i.split('_')[2]), int(i.split('_')[1][1:])) for i in data.columns],
        names=['reps', 'time'])

    data.columns = idx

    return data, names

def _cov(x,y, axis=0):
    return np.ma.mean(x*y, axis=axis) - (np.ma.mean(x, axis=axis)*np.ma.mean(y, axis=axis))

def prep_input(abundance_file, counts_file):
    ab = _load_abundance_thresholds(abundance_file)
    tps = _load_timepoint_counts(counts_file)

    _tps, names = _convert_timepoint_counts(tps)
    _abundance = _convert_abundance_thresholds(ab)

    names.loc[names['target_a_id'].str.contains('NonTargeting'), 'target_a_id'] = '0'
    names.loc[names['target_b_id'].str.contains('NonTargeting'), 'target_b_id'] = '0'

    good = ~(names['target_a_id'] == names['target_b_id'])
    good_data = _tps.loc[good]
    good_names = names.loc[good]

    cpA = good_names['probe_a_id'].apply(lambda x: '0' + x if 'NonTargeting' in x else x)
    cpB = good_names['probe_b_id'].apply(lambda x: '0' + x if 'NonTargeting' in x else x)

    pswitch = cpA > cpB
    phold = cpA.loc[pswitch]
    cpA.loc[pswitch] = cpB.loc[pswitch]
    cpB.loc[pswitch] = phold

    probes = pd.concat([cpA, cpB]).unique()
    probes.sort()
    nprobes = len(probes)

    cgA = good_names['target_a_id']
    cgB = good_names['target_b_id']

    genes = pd.concat([cgA, cgB]).unique()
    genes.sort()

    n = genes.shape[0]
    mm = n*(n-1)/2

    gswitch = cgA > cgB
    ghold = cgA.loc[gswitch]

    cgA_c = cgA.copy() # Avoid the copy warning
    cgB_c = cgB.copy()

    cgA_c.loc[gswitch] = cgB.loc[gswitch]
    cgB_c.loc[gswitch] = ghold

    cgA = cgA_c
    cgB = cgB_c

    gA_gB = cgA.str.cat(cgB, sep='_')
    pA_pB = cpA.str.cat(cpB, sep='_')

    good_data.values[good_data.values == 0] = 1
    abundance = good_data.sum(axis=0)
    y = np.log2(good_data/abundance)

    ab0 = pd.Series(ab.values.ravel() - np.log2(abundance.values), index=abundance.index)

    counts = np.array([y[1].values, y[2].values])
    ab = np.array([ab0[1].values, ab0[2].values])[...,np.newaxis].transpose(0,2,1)

    return ab, counts, good_names

def _validate_time(counts_shape, times):
    n_reps, n_samples, n_timepts = counts_shape

    if n_timepts not in times.shape:
        raise ValueError('Inconsistent times dimensions!')

    if len(times.shape) > 3:
        raise ValueError('Times have too many dimensions')

    elif len(times.shape) == 1:
        times = repmat(times, n_samples, n_reps).reshape(counts_shape)

    elif len(times.shape) == 2:
        times = repmat(times, n_samples, 1).reshape(counts_shape)

    elif len(times.shape) == 3:
        times = times.rehsape(counts_shape)

    return times

# def _validate_abundance(counts_shape, abundance):
#     n_reps, n_samples, n_timepts = counts_shape
#     a,b = abundance.shape
#
#     if a not in counts_shape or b not in counts_shape:
#         raise ValueError('Inconsistent abundance dimensions!')
#
#     return abundance.reshape((n_reps, n_timepts))


def fit_ac_fc(abundance_file, counts_file, times, n_good=2):
    '''fit_ac_fc

    This is a line by line recapitulation of Amanda's code (vectorized using
    matrix multiplication implemented in python).

    TODO: Figure out why lfdr is different
    TODO: Reimplement lfdr in python
    TODO: Refactor
    TODO: Confirm input output formats

    Note: validation is based on the counts_file'''

    ab, counts, names = prep_input(abundance_file, counts_file)
    n_reps, n_samples, n_timepts = counts.shape

    times = _validate_time(counts.shape, times)
    #ab = _validate_abundance(counts.shape, ab)

    bad = (counts > ab).sum(axis=2) < n_good
    allbad = bad.all(axis=0)

    mask = counts > ab
    useless = mask.sum(axis=2) < 2 #NOTE: Not sure this is needed
    mask[useless] = False

    counts_masked = np.ma.array(data=counts, mask=~mask)
    time_masked = np.ma.array(data=times, mask=~mask)

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

    df = mask.sum(axis=2).sum(axis=0) - 2
    df[allbad] = 0

    counts_masked_2d = np.ma.array(data=np.hstack(counts), mask=~np.hstack(mask))
    xfit_masked_2d = np.ma.array(data=np.hstack(xfit), mask=~np.hstack(mask))
    time_masked_2d = np.ma.array(data=np.hstack(times), mask=~np.hstack(mask))

    n = xfit_masked_2d - counts_masked_2d
    num = np.sqrt(np.power(n,2).sum(axis=1))

    d = time_masked_2d - time_masked_2d.mean(axis=1)[...,np.newaxis]
    denom = np.sqrt(np.power(d, 2).sum(axis=1))

    #NOTE: Hardcode
    sdfc = np.divide(num, denom).data
    sdfc[allbad] = 0.1

    has_sd = df > 0
    median_sd = np.median(sdfc[has_sd])
    sdfc[~has_sd] = median_sd

    df_masked = np.ma.array(data=df, mask=~(df > 0))

    tstat_masked = np.ma.divide(fc, np.ma.divide(sdfc, np.ma.sqrt(df_masked)))
    tstat_masked.data[tstat_masked.mask] = 1

    tstat = tstat_masked.data
    p_t = np.array([2*t.cdf(-abs(i),j) if j > 0 else 1 for i,j in zip(tstat, df)])

    qvalue = importr('qvalue')
    r_pt = numpy2ri(p_t[has_sd]) #Missing has_sd

    lfdr = np.ones((n_samples))
    lfdr[has_sd] = np.array(qvalue.lfdr(r_pt, method="bootstrap")) #NOTE: Still deviates from ctg

    #return ac, fc, allbad, sdfc, df, p_t, lfdr, names, lmbda, xfit, mask
    return ac, fc, allbad, sdfc, df, p_t, lfdr, names


if __name__ == "__main__":
    abundance_file = os.path.join(config.A549_test, "A549_abundance_thresholds.txt")
    counts_file = os.path.join(config.A549_test, "A549_timepoint_counts.txt")
    times = np.array([[3,14, 21, 28], [3,14,21,28]])

    fit_ac_fc(abundance_file, counts_file, times)
