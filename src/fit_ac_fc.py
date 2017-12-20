"""Initializes abundance and fitnesses for the iterative least squares

This module initializes abundance and fitness for iterative least squares later.

TODO:
    * Update docstrings
    * Separate conversion tools (to change existing data structure to new ones)
    * Change variable names to something more descriptive

"""

import os
import pandas as pd

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

if __name__ == "__main__":
    from config import config

    ab = load_abundance_thresholds(os.path.join(config.A549_test, \
            "A549_abundance_thresholds.txt"))
    tps = load_timepoint_counts(os.path.join(config.A549_test, \
            "A549_timepoint_counts.txt"))

    #print(convert_abundance_thresholds(ab))

    print(convert_timepoint_counts(tps))
