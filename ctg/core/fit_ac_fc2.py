import os
import sys
from collections import namedtuple

import numpy as np
from numpy.matlib import repmat
import pandas as pd
from scipy.stats import t

from ctg.core.config import config
import ctg.core.calculate_abundance as calculate_abundance

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
        min_good_tpts=2,
    ): 
        
        self.data = dataframe
        self.mask = mask
        self.min_good_tpts = min_good_tpts

        #Separating names and timepoint data
        if names is None: 
            self.names = dataframe.iloc[:, :col]
            self.data = dataframe.iloc[:, col:]
        else: 
            self.names = names
        
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

        #Log Transforming counts
        good_data = self.data.loc[good]
        good_data.values[good_data.values == 0] = 1
        abundance = good_data.sum(axis=0)
        y = np.log2(good_data/abundance)

        self.data = y
        self.good_names = good_names

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

            if len(arr) != 3 or arr[1][0] != "T": 
                raise ValueError("Column headers are expected to be in the form {NAME}_T{DAYS}_{REP}")
            
            container.append((int(arr[2]), int(arr[1][1:]))) # (Rep, Timepoint)

        reps = set([i[0] for i in container])
        if reps != set(range(1, len(reps) + 1)): 
            raise ValueError("Expect reps to be integers starting from 1.")

        self.n_reps = len(reps) 
        self.timepoints = [[] for _ in range(self.n_reps)]
        indexes = [[] for _ in range(self.n_reps)]

        for ind, tup in enumerate(container): 
            i,j = tup
            indexes[i - 1].append(ind)
            self.timepoints[i - 1].append(j)

        self.data_indexes = indexes


    def add_mask(self): 
        """Creates a mask for which timepoints did not meet abundance threshold""" 

        if not hasattr(self, "abundance_threshold"): 
            raise ValueError("Cannot create mask without abundance threshold!")

        Masker = namedtuple("Masker", ["mask", "bad", "allbad"])

        mask = self.data.values > self.abundance_threshold.values

        bad = [mask[:, index].sum(axis=1) < self.min_good_tpts \
            for index in self.data_indexes
        ]

        bad = np.vstack(bad)
        allbad = bad.all(axis=0)

        for i, index in enumerate(self.data_indexes): 
            mask[:, index][bad[i]] = False

        self.mask = Masker(mask, bad, allbad)


    def add_abundance_threshold(
        self, 
        abundance_thresholds=None, 
        min_counts_threshold=10
    ): 
        """Calculate or add abundance threshold"""

        if abundance_thresholds is None: 
            self.abundance_thresholds = calculate_abundance.calculate_abundance(
                self.data.values, 
                min_counts_threshold=min_counts_threshold
            )

        else: 
            if set(self.data.columns.tolist()) != set(abundance_thresholds.index.tolist()): 
                raise ValueError("Expected abundance threhsolds to have index the same as the data column.")
            
            self.abundance_threshold = abundance_thresholds.T

        return self 


    def calculate_construct_fitness(self): 
        """Calculate construct fitness"""

        fitnesses = []
        var_times = []
        times_list = []

        mean_time_list = []
        mean_counts_list = []

        for index, times in zip(self.data_indexes, self.timepoints):     
            counts = self.data.iloc[:, index].values
            mask = ~self.mask.mask[:, index]

            times = repmat(times, counts.shape[0], 1)
            times_list.append(times)

            counts_masked = np.ma.array(data=counts, mask=mask)
            time_masked = np.ma.array(data=times, mask=mask)

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
        #TODO: mask the divide so you can supress the warning
        self.fitness = np.divide(np.sum(fitnesses, axis=0), np.sum(var_times, axis=0))
        self.fitness[self.mask.allbad] = 0

        self.f = fitnesses
        self.var_times = var_times 

        acs = []
        xfits = []

        rep = 0
        for mean_counts, mean_time in zip(mean_time_list, mean_counts_list): 
            ac = mean_counts - (self.fitness*mean_time)
            ac[self.mask.bad[rep]] = counts[self.mask.bad[rep], 0] #Take the first count

            alpha = -np.log2(np.power(2,ac).sum())
            ac = ac + alpha
        
            lmbda = -np.log2(
                np.power(
                    2, 
                    ac[..., np.newaxis] + self.fitness[..., np.newaxis]*times[rep]
                ).sum()
            )
            xfit = ac[...,np.newaxis] + self.fitness[...,np.newaxis]*times[rep] + lmbda
        
            acs.append(ac)
            xfits.append(xfit)

            rep += 1

        self.acs = acs
        self.xfits = xfits
        self.times_list = times_list

        return self


    def  estimate_pvalue(self): 
        """Estimate p-value"""

        if not hasattr(self, "fitness"): 
            raise RuntimeError("Please calculate construct fitness first using calculate_construct_fitness")
    
        mask = ~self.mask.mask

        counts_masked_2d = np.ma.array(data=self.data.values, mask=mask)
        xfit_masked_2d = np.ma.array(data=np.hstack(self.xfits), mask=mask)
        time_masked_2d = np.ma.array(data=np.hstack(self.times_list), mask=mask)
        
        #Degrees of freedom
        #df = mask.mask.sum(axis=2).sum(axis=0) - 2 #Hard-coded (2?)
        df = self.mask.mask.sum(axis=1) - 2
        df[self.mask.allbad] = 0

        n = xfit_masked_2d - counts_masked_2d
        num = np.sqrt(np.power(n,2).sum(axis=1))

        d = time_masked_2d - time_masked_2d.mean(axis=1)[...,np.newaxis]
        denom = np.sqrt(np.power(d, 2).sum(axis=1))

        #NOTE: Hardcode
        sdfc = np.divide(num, denom).data
        sdfc[self.mask.allbad] = 0.1

        has_sd = df > 0
        median_sd = np.median(sdfc[has_sd])
        sdfc[~has_sd] = median_sd

        df_masked = np.ma.array(data=df, mask=~(df > 0))

        tstat_masked = np.ma.divide(self.fitness, np.ma.divide(sdfc, np.ma.sqrt(df_masked)))
        tstat_masked.data[tstat_masked.mask] = 1

        tstat = tstat_masked.data
        p_t = np.array([2*t.cdf(-abs(i),j) if j > 0 else 1 for i,j in zip(tstat, df)])
        
        self.sdfc = sdfc 
        self.p_t = p_t 

        return self 

    
    def fit_ac_fc(self, abundance_df): 
        """Wrapper method to run the entire pipeline""" 
        
        self.add_abundance_threshold(abundance_df)
        self._sanitize_names()
        self.add_mask() 
        self.calculate_construct_fitness()
        self.estimate_pvalue()

        return self.acs, self.fitness, self.mask.allbad, self.sdfc, self.p_t, self.good_names