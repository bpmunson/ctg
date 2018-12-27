import os
import sys

import numpy as np
from numpy.matlib import repmat
import pandas as pd
from scipy.stats import t

from ctg.core.config import config
#import config
import ctg.core.calculate_abundance as calculate_abundance
#import calculate_abundance

# import warnings
# from pandas.core.common import SettingWithCopyWarning
# warnings.simplefilter('error', SettingWithCopyWarning)

'''
The file format below is hardcoded for Amanda's files. This can be changed by
later by passing a global config object for the run (or just additional
function argumments).

Currently, the time points should be of the format

construct_id    probe_a_id      probe_b_id      target_a_id     target_b_id     {NAME}_T{DAYS}_{REP}

In addition, reps are defined to be the first set of levels in the loading functions

'''

class masker(object):
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


def ma_cov(x,y, axis=0):
    """Calculates the covariance from masked numpy arrays"""

    return np.ma.mean(x*y, axis=axis) - (np.ma.mean(x, axis=axis)*np.ma.mean(y, axis=axis))

class Counts(object): 
    """Counts object to contain information related to the counts"""

    def __init__(
        self, 
        dataframe, 
        mask=None,
        names=None, 
        col=5,
        min_good_tpts=0,
    ): 
        
        self.data = dataframe
        self.mask = mask
        self.min_good_tpts = min_good_tpts

        #Separating names and timepoint data
        if names is not None: 
            self.names = dataframe.iloc[:, :col]
            self.data = dataframe.iloc[:, col:]
        else: 
            self.names = names
        
        self._sanitize_names()
        self._parse_header()


    @classmethod
    def from_file(
        cls, 
        file, 
        names=None, 
        col=5, 
        sep='\s+',
        **kwargs
    ): 
        """Constructing an Counts object using a text file"""

        kwargs['sep'] = sep # Setting default
        df = pd.read_csv(file, **kwargs)

        return Counts(
            df,  
            names=names, 
            col=col,
        )

    def _sanitize_names(self): 
        """Leftover compliance from previous pipeline"""
        if self.names is None: 
            raise RuntimeError("Cannot sanitize without providing names!") 

        good = ~(self.names['target_a_id'] == self.names['target_b_id'])
        good_names = self.names.loc[good]

        """
        #I think chunk of code is useless 
        
        cpA = good_names['probe_a_id']
        cpB = good_names['probe_b_id']

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

        gswitch = cgA > cgB
        ghold = cgA.loc[gswitch]

        # had to remove the copy to ensure probes and target names matches
        cgA_c = cgA
        cgB_c = cgB 


        cgA_c.loc[gswitch] = cgB.loc[gswitch]
        cgB_c.loc[gswitch] = ghold

        cgA = cgA_c
        cgB = cgB_c

        gA_gB = cgA.str.cat(cgB, sep='_')
        pA_pB = cpA.str.cat(cpB, sep='_')
        """

        #Log Transforming counts
        good_data = self.data.loc[good]
        good_data.values[good_data.values == 0] = 1
        abundance = good_data.sum(axis=0)
        y = np.log2(good_data/abundance)

        self.data = y

        if hasattr(self, "abundance_threshold"): 
            self.abundance_threshold = pd.Series(
                self.abundance_threshold.values.ravel() - np.log2(abundance.values), 
                index=abundance.index
            )

    def _parse_header(self): 
        """Extract number of replicates and timepoints from header"""
        container = []
        for i in self.data.columns: 
            arr = i.split('_') 

            if len(arr) != 2 or arr[1][0] != "T": 
                raise ValueError("Column headers are expected to be in the form {NAME}_T{DAYS}_{REP}")
            
            container.append((arr[2], int(arr[1][1:]))) # (Rep, Timepoint)

        reps = set([i[0] for i in container])
        if reps != set(range(1, len(reps) + 1)): 
            raise ValueError("Expect reps to be integers starting from 1.")

        self.n_reps = len(reps) 
        self.timepoints = [[] for _ in self.n_reps]
        indexes = [[] for _ in self.n_reps]

        for ind, tup in enumerate(container): 
            i,j = tup
            indexes[i - 1] = ind
            self.timepoints[i - 1] = j

        self.data_indexes = indexes


    def add_abundance_threshold(
        self, 
        abundance_thresholds=None, 
        min_counts_threshold=10
    ): 
        """Calculate or add abundance threshold"""

        if abundance_thresholds is None: 
            self.abundance_thresholds = calculate_abundance.calculate_abundance(
                _tps, 
                min_counts_threshold=min_counts_threshold
            )

        else: 
            if set(self.data.columns.tolist()) != (abundance_thresholds.index.tolist()): 
                raise ValueError("Expected abundance threhsolds to have index the same as the data column.")
            
            self.abundance_threshold = abundance_thresholds

        # Add abundance threshold as a mask 
        self.mask = masker(
            self.data.values, 
            self.abundance_threshold, 
            self.min_good_tpts
        )

        return self


    def calculate_construct_fitness(self): 
        """Calculate construct fitness"""

        fitnesses = []
        var_times = []
        mean_time_list = []
        mean_counts_list = []
        for index, times in zip(self.data_indexes, self.timepoints):     
            counts = self.data.iloc[:, index].values

            counts_masked = np.ma.array(data=counts, mask=~self.mask.mask)
            time_masked = np.ma.array(data=times, mask=~self.mask.mask)

            #Collapsing across timepoints
            mean_counts = np.ma.mean(counts_masked, axis=1).data
            mean_counts_list.append(mean_counts)
            mean_time = np.ma.mean(time_masked, axis=1).data
            mean_time_list.append(mean_time)

            var_time = np.ma.var(time_masked,axis=1).data
            f = ma_cov(counts_masked, time_masked, axis=1).data

            var_times.append(var_time)
            fitnesses.append(f)

        #Collapsing across replicates 
        self.fitness = np.divide(np.sum(fitnesses), np.sum(var_times))
        self.fitness[self.mask.allbad] = 0

        #TODO:THIS PART NEEDS ATTENTION. INCONSISTENT DEFINITION OF AC (There
        #should be two of them?)

        acs = []
        xfits = []
        for mean_counts, mean_time in zip(mean_time_list, mean_counts_list): 
            ac = mean_counts - (self.fitness*mean_time)
            ac[mask.bad] = counts[mask.bad,0]

            alpha = -np.log2(np.power(2,ac).sum(axis=1))[...,np.newaxis]
            ac = ac + alpha
        
            lmbda = -np.log2(np.power(2, ac[...,np.newaxis] + fc[np.newaxis,...,np.newaxis]*times).sum(axis=1))
            xfit = ac[...,np.newaxis] + fc[np.newaxis,...,np.newaxis]*times + lmbda[...,np.newaxis]
        
            xfits.append(xfit)

        return self


    def  estimate_pvalue(self): 
        """Estimate p-value"""

        if not hasattr(self, "construct_fitness"): 
            raise RuntimeError("Please calculate construct fitness first using calculate_construct_fitness")
    
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

    
    def fit_ac_fc(self): 
        self.add_abundance_threshold()
        self.calculate_construct_fitness()
        self.estimate_pvalue()