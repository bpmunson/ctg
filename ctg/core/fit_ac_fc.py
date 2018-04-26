import os
import sys

import numpy as np
from numpy.matlib import repmat
import pandas as pd
from scipy.stats import t

from ctg.core.config import config
import ctg.core.calculate_abundance as calculate_abundance


'''
The file format below is hardcoded for Amanda's files. This can be changed by
later by passing a global config object for the run (or just additional
function argumments).

Currently, the time points should be of the format

construct_id    probe_a_id      probe_b_id      target_a_id     target_b_id     {NAME}_T{DAYS}_{REP}

In addition, reps are defined to be the first set of levels in the loading functions

'''

def _parse_index(idx): 
    names = []
    for i in idx: 
        arr = i.split('_')
        names.append((int(arr[2]), int(arr[1][1:])))

    return names

def _reshape_using_multiindex(ary, idx): 
    if sorted(list(idx)) != list(idx): 
        raise ValueError("idx needs to be a sorted list")
        
    unique_idx = []
    unique = []
    for i,j in enumerate(idx): 
        if j not in unique: 
            unique.append(j)
            unique_idx.append(i)
            
    return np.stack(np.split(ary, unique_idx[1:], axis=1))

def _read_abundance_thresholds_file(fn, sep='\s+', **kwargs):
    """Loads abundance thresholds.

    """

    ab = pd.read_csv(fn, sep=sep, index_col=0, **kwargs)
    names = _parse_index(ab.index)

    idx = pd.MultiIndex.from_tuples(names, names=['reps', 'time'])
    ab.index = idx
    ab.sort_index(axis=0, inplace=True)

    return ab

def _read_counts_file(fileName, n_reps=None, t=None, columns_map=None, method='implicit', col=5, sep='\s+',**kwargs): 
    ''' Reads counts file

    All other kwargs are passed onto pandas read_csv
    sep will be added to the kwargs so it will also be passed into read_csv
    
    col is the index column of the start of the counts
    '''
    
    kwargs['sep'] = sep # Make sure sep is also passed

    #Input validation
    if method not in ['explicit', 'implicit']: 
        raise ValueError('method must only be explicit or implicit')

    elif method == 'explicit' and (n_reps is None or t is None or columns_map is None): 
        raise ValueError('n_reps (number of replicates), t (time points), and columns_map must be supplied if method is explicit!')

    #Read file as dataframe
    df = pd.read_csv(fileName, **kwargs)
    prefix = df.iloc[:, :col]
    counts = df.iloc[:, col:]
    
    #Decide the column mappings
    if method == 'explicit': 
        #Validating optional arguments
        expected_columns = n_reps*len(t)
        if expected_columns != counts.shape[1]: 
            raise ValueError("Expected %s columns. Got %s" % (expected_columns, counts.shape[1]))
            
        if len(columns_map) != len(t): 
            raise ValueError("Each time point needs a list of column indices!")
            
        #Creating tuples for multiindex
        names = [None for _ in range(n_reps*len(t))]
        for i_ind, i in enumerate(columns_map): 
            if len(i) != n_reps: 
                raise ValueError("Each list within columns_map needs to have the same length as n_reps")
                
            for j_ind, j in enumerate(i): 
                names[j] = (j_ind + 1, t[i_ind])
            
    elif method == 'implicit':         
        names = _parse_index(counts.columns)
            
    #Reshape the counts dataframe
    idx = pd.MultiIndex.from_tuples(names, names=['reps', 'time'])
    counts.columns = idx    
    counts.sort_index(axis=1, inplace=True)
    
    counts_mat = _reshape_using_multiindex(counts.values, counts.columns.labels[0])

    if method == 'explicit': 
        try: 
            assert counts_mat.shape[0] == n_reps
            assert counts_mat.shape[2] == len(t)
        
        except AssertionError: 
            print(counts_mat.shape)
            raise RuntimeError("Internal error... Counts matrix was reshaped incorrectly!")
    
    return prefix, counts, counts_mat

def _cov(x,y, axis=0):
    return np.ma.mean(x*y, axis=axis) - (np.ma.mean(x, axis=axis)*np.ma.mean(y, axis=axis))

def _validate_counts(counts, replicate_axis=0, samples_axis=1, timepoints_axis=2):
    if not isinstance(counts, np.ndarray):
        raise ValueError('Please input a numpy array for counts!')

    if counts.ndim not in [2,3]:
        raise ValueError('counts need to be 2D or 3D numpy array!')

    if counts.ndim == 2:
        counts = counts.reshape((1, counts.shape[samples_axis], counts.shape[timepoints_axis]))

    elif counts.ndim == 3:
        s = counts.shape
        new_shape = (s[replicate_axis], s[samples_axis], s[timepoints_axis])
        counts = counts.reshape(new_shape)

    return counts

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

class _masker(object):
    '''Object to generate and contain the mask needed to separate the good and bad counts.

    Args: 
        counts (numpy array): counts for each constructs
        abundance (numpy array): abundance threshold for counts (needs to be broadcastable with counts)
        min_good_tpts (int, optional): minimum of number of good timepoints to be considered "good"

    '''
    def __init__(self, counts, abundance, min_good_tpts=2): 
        self.mask = counts > abundance
        self.bad = self.mask.sum(axis=2) < min_good_tpts
        self.allbad = self.bad.all(axis=0)

        self.mask[self.bad] = False

def _prep_input(abundance_file, counts_file, names=None, n_reps=None, t=None, columns_map=None, method='implicit', col=5, sep='\s+', min_counts_threshold=10, **kwargs):
    if isinstance(counts_file, str):
        names, _tps, _tps_mat = _read_counts_file(counts_file, n_reps=n_reps, t=t, columns_map=columns_map, method=method, col=col, sep=sep, **kwargs)

    elif isinstance(counts_file, pd.DataFrame):
        _tps = counts_file #TODO: This option needs work (ensure multi-index, build?)

    if isinstance(abundance_file, str):
        ab = _read_abundance_thresholds_file(abundance_file) 

    elif abundance_file is None:
        #TODO: Yes this is a hack...
        ab = calculate_abundance.calculate_abundance(_tps, min_counts_threshold=min_counts_threshold)

    elif isinstance(abundance_file, pd.DataFrame) or isinstance(abundance_file, pd.Series):
        #TODO: Verify that series work also

        ab = abundance_file

    else:
        raise ValueError('Expected abundance to be of type {str, pd.DataFrame, pd.Series, None}. Received %s instead' % type(abundance_file))

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

    counts = np.array([y[i].values for i in y.columns.levels[0]]) # Assume 'reps' to be the first set of levels (see above)
    ab = np.array([ab0[i].values for i in ab0.index.levels[0]])[...,np.newaxis].transpose(0,2,1) #ditto

    return ab, counts, good_names    

def _estimate_fitness(counts, times, mask): 
    '''Estimate the fitness based on counts and time.

    Performs a regression by estimating the counts covariance (across time) and the variance time. 
    
    Args: 
        counts (numpy array): Counts for each construct
        times (numpy array): Time for each column of counts
        mask (_masker object): contains the boolean mask of counts
        
    Returns: 
        fc (numpy array): fitness estimate for each construct
        ac (numpy array): initial count for each construct
        xfit (numpy array): estimated fitness based on the regression for each construct
    '''
    counts_masked = np.ma.array(data=counts, mask=~mask.mask)
    time_masked = np.ma.array(data=times, mask=~mask.mask)

    mean_counts = np.ma.mean(counts_masked, axis=2).data
    mean_time = np.ma.mean(time_masked, axis=2).data
    var_time = np.ma.var(time_masked,axis=2).data
    f = _cov(counts_masked, time_masked, axis=2).data

    fc = np.divide(f.sum(axis=0), var_time.sum(axis=0))
    fc[mask.allbad] = 0

    ac = mean_counts - (fc*mean_time)
    ac[mask.bad] = counts[mask.bad,0]

    alpha = -np.log2(np.power(2,ac).sum(axis=1))[...,np.newaxis]
    ac = ac + alpha
    
    lmbda = -np.log2(np.power(2, ac[...,np.newaxis] + fc[np.newaxis,...,np.newaxis]*times).sum(axis=1))
    xfit = ac[...,np.newaxis] + fc[np.newaxis,...,np.newaxis]*times + lmbda[...,np.newaxis].transpose(0,2,1)
    
    return fc, ac, xfit

def _estimate_pvalue(xfit, fc, counts, times, mask):    
    '''Estimates the confidence in each regression

    Estimates a p-value for each construct denoting the confidence for the estimated fitness.

    Args: 
        xfit (numpy array): estimated fitness based on the regression for each construct
        fc (numpy array): fitness estimate for each construct
        counts (numpy array): Counts for each construct
        times (numpy array): Time for each column of counts
        mask (_masker object): contains the boolean mask of counts

    Returns: 
        sdfc (numpy array): standard deviations of the fitness
        p_t (numpy array): p-value associated with each construct
    '''
    counts_masked_2d = np.ma.array(data=np.hstack(counts), mask=~np.hstack(mask.mask))
    xfit_masked_2d = np.ma.array(data=np.hstack(xfit), mask=~np.hstack(mask.mask))
    time_masked_2d = np.ma.array(data=np.hstack(times), mask=~np.hstack(mask.mask))
    
    df = mask.mask.sum(axis=2).sum(axis=0) - 2
    df[mask.allbad] = 0

    n = xfit_masked_2d - counts_masked_2d
    num = np.sqrt(np.power(n,2).sum(axis=1))

    d = time_masked_2d - time_masked_2d.mean(axis=1)[...,np.newaxis]
    denom = np.sqrt(np.power(d, 2).sum(axis=1))

    #NOTE: Hardcode
    sdfc = np.divide(num, denom).data
    sdfc[mask.allbad] = 0.1

    has_sd = df > 0
    median_sd = np.median(sdfc[has_sd])
    sdfc[~has_sd] = median_sd

    df_masked = np.ma.array(data=df, mask=~(df > 0))

    tstat_masked = np.ma.divide(fc, np.ma.divide(sdfc, np.ma.sqrt(df_masked)))
    tstat_masked.data[tstat_masked.mask] = 1

    tstat = tstat_masked.data
    p_t = np.array([2*t.cdf(-abs(i),j) if j > 0 else 1 for i,j in zip(tstat, df)])
    
    return sdfc, p_t

def fit_ac_fc(abundance, counts, times, method='implicit', n_reps=None, columns_map=None, col=5, 
                names=None, min_good_tpts=2, min_counts_threshold=10,):
                #replicate_axis=0, samples_axis=1, timepoints_axis=2,
                #verbose=False): 
    
    '''fit_ac_fc

    Fits a fitness (growth rate) to the counts vs timepoints. After normalization, the function simply does
    a linear regression to to across the timepoints (while masking away the "bad" timepoints, i.e. those that
    do not have enough counts).

    Args: 
        abundance (str, DataFrame, Series, None): The abundance threshold to be considered "good".
        counts (str, DataFrame): The counts for each construct
            Can be string (file location) or pre-parsed DataFrame. If it is a DataFrame, it needs to 
            have Multi-Index as column names (first level is the replicates (integer 1,2,...), and 
            second level is the time points (integer 3,7,...))
        times (numpy array): Time points for the counts vs time regression
        method (str): either implicit or explicit. If implicit, the column names of the counts file is read, 
            and the number of replicates and time points are inferred. Otherwise, the counts dataframe will be 
            forced into the shape according to the given parameters.
        n_reps (int, optional): Number of replicates. This argument is required if the method is explicit.
        columns_map (list of list of int, optional): Provides the explicit column to time points and replicates 
            mapping. The list should have the same length has the number of time points. Each element contains
            a list whose elements are column index (0 index) in the order of the replicate number.
        col (int, optional): The column index where the counts start. Any columns before this is considered
            the names dataframe
        names (dataframe, optional): If counts were passed as a DataFrame, names must also be passed. This 
            DataFrame contains the construct id, probe, id, etc (see top doctsring)
        min_good_tpts (int, optional): Minimum number of timepoints that pass abundance threshold to be 
            considered good. 
        min_counts_threshold (int, optional): Minimum number of counts required to be considered good
        replicate_axis

    Returns: 
        ac (numpy array): initial count for each construct
        fc (numpy array): fitness for each construct
        mask.allbad (numpy array): boolean array for each construct. True if the construct is deemed bad
        sdfc (numpy array): standard deviations of the fitness estimate for each construct
        p_t (numpy array): p-value for each construct 
        names (dataframe): All the prefix information (see top doctsring)      
    '''
    
    abundance, counts, names = _prep_input(abundance, counts,names=names,
                                        t=times, n_reps=n_reps, columns_map=columns_map, method=method, col=col,
                                        min_counts_threshold=min_counts_threshold,
                                        verbose=verbose)

    counts = _validate_counts(counts)

    n_reps, n_samples, n_timepts = counts.shape
    ab = abundance #TODO: Change ab

    times = _validate_time(counts.shape, times)
    #ab = _validate_abundance(counts.shape, ab)

    mask = _masker(counts, ab, min_good_tpts=min_good_tpts)
    fc, ac, xfit = _estimate_fitness(counts, times, mask)
    sdfc, p_t = _estimate_pvalue(xfit, fc, counts, times, mask)
    
    return ac, fc, mask.allbad, sdfc, p_t, names

if __name__ == "__main__":
    # abundance_file = os.path.join(config.A549_test, "A549_abundance_thresholds.txt")
    # counts_file = os.path.join(config.A549_test, "A549_timepoint_counts.txt")
    # times = np.array([[3,14, 21, 28], [3,14,21,28]])
    #
    # fit_ac_fc(abundance_file, counts_file, times)


    abundance_file = os.path.join(config.A549_test, "A549_abundance_thresholds_rep1.txt")
    counts_file = os.path.join(config.A549_test, "A549_timepoint_counts_rep1.csv")
    times = np.array([3,14,21,28])

    #prep_input(abundance_file, counts_file)
    fit_ac_fc(abundance_file, counts_file, times)
