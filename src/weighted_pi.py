"""
Docstring


This is a rewrite of the pi score weighting from Roman/Amanda's CTG pacakge.
Ported to pure python.

Todo:
    method documentation
    think about alt weighting strategies
    unit tests

    get rid of that regex garbage
    rethink pandas strategy
"""
import re
import os
import argparse
import pandas as pd
import numpy as np



def _get_null_dataframe(df, indices=None, null_regex=".*NonTargeting.*"):
    """ Get the null dataframe corresponding to non targeting probes
    """
    # if no index range was included use regex 
    if indices is None:
        r = re.compile(null_regex)
        probes = [i for i, probe in enumerate(df.index) if r.match(probe)]

    # calculate mean of fitnesses for null probes
    null_df = df.iloc[probes, :]

    return null_df

def _get_null_indices(df, indices=False, null_regex=".*NonTargeting.*"):
    """ Get the null indices of the dataframe corresponding to non targeting probes
    Args:
        df (dataframe): Pandas dataframe to get the probe names or indices from, row index should contains probe names
        indices (bool): True to return indices of null probes in df, False to return probe names of null probes in df
        regex (str): a python regular expression to match null probes; default dual crispr style is '.*NonTargeting.*'
    Returns:
        ids (list): a list of probe names corresponding to the null target
    """

    r = re.compile(null_regex)
    ids = [i for i, probe in enumerate(df.index) if r.match(probe)]

    return ids

def _get_null_probes_names(df, null_regex=".*NonTargeting.*"):
    """ Get the null probe names of the dataframe corresponding to non targeting probes
    Args:
        df (dataframe): Pandas dataframe to get the probe names or indices from, row index should contains probe names
        regex (str): a python regular expression to match null probes; default dual crispr style is '.*NonTargeting.*'
    Returns:
        ids (list): a list of probe names corresponding to the null target
    """

    r = re.compile(null_regex)
    ids = [probe for i, probe in enumerate(df.index) if r.match(probe)]
    
    return ids

def _build_target_to_probe_table(fn, sep="\t", header=0, comment = "#", null_regex=".*NonTargeting.*"):
    """ Loads library definition
        Todo: this assumes a CTG style library definition, should probably abstract it
    """
    # read in library definitions file
    library = pd.read_csv(fn, sep=sep, header=header, comment=comment)
    # take only the target and probe id columns
    a = library[['target_a_id','probe_a_id']]
    b = library[['target_b_id', 'probe_b_id']]
    # relabel with identical headers
    a.columns = ['target_id', 'probe_id']
    b.columns = ['target_id', 'probe_id']
    # concatenate together
    c = pd.concat([a,b])
    # uniquify
    c = c.drop_duplicates()
    c.index = c['probe_id']
    c.drop('probe_id', axis=1, inplace=True)

    # replace NonTargetingNNNN target id 
    null_probes = _get_null_probes_names(c, null_regex = null_regex )
    c.loc[null_probes, 'target_id'] = "0NonTargeting" # TODO: change this hardcoding
    return c

def rank_probes(fp, library, null=True, null_probes = None, null_regex = ".*NonTargeting.*"):
    """Rank the probes by 
    Todo:
    should we require the library_definition?
    add special case for null probes where we want the closest ones to zero
    """

    # add target ids to the fitness dataframe
    fpc = pd.merge(fp, library, left_index = True, right_index=True)



    # get absolute value of fitness
    fpc.loc[:,'fitness_abs'] = fpc.loc[:,'fitness'].abs()

    # get rank
    # for the knockout probesthe assumption here is that we want the probes which deviate furthest from zero 
    fpc.loc[:,'rank'] = fpc.groupby('target_id')['fitness_abs'].rank(method='max', ascending=False)

    # do something special for the null probes
    if null:
        if null_probes is None:
            # get the null probe names
            null_probes = _get_null_probes_names(fpc, null_regex = null_regex )

        # replace the ranks for the null probes
        # for the null probes we want to take the probes which are closest to zero 
        fpc.loc[null_probes,'rank'] = fpc.loc[null_probes, :].groupby('target_id')['fitness_abs'].rank(method='min', ascending=True)
 

    return fpc

def ansatz_construct_weights(fpr, n_probes_per_target=2):
    """Build the construct weight matrix given a dataframe of probe ranks by target

    Args: 
        fpr (dataframe): a pandas dataframe with indices of probe names and a column corresponding to 
                        single probe ranks with column header 'rank'. may contain extra columns.
    Returns:
        construct_weights (dataframe): a NxN dataframe of rank based weights for each probe pair.

    """

    # get number of probes

    n = fpr.shape[0]
    construct_weights = np.zeros((n,n))
    rank_iloc = list(fpr.columns).index('rank')
    for i in range(n):
        for j in range(n):
            # so
            # rank(1) * rank(1) is the best and gets a weight of 4
            # rank(2) * rank(1) is second best and gets a weight of 2
            # rank(2) * rank(2) is next and gets a weight of 1
            # rank(3) * rank(3) is worst and we wont use it so it gets a rank of 0

            # abstract this for any number of probes per construct
            construct_weights[i][j] =   (fpr.iloc[i,rank_iloc] - (n_probes_per_target+1)) * \
                                        (fpr.iloc[j,rank_iloc] - (n_probes_per_target+1))

    construct_weights = pd.DataFrame(construct_weights, index=fpr.index, columns=fpr.index)
    return construct_weights

def ansatz_target_fitness(fpr, n_probes_per_target=2):
    """Build the construct weight matrix given a dataframe of probe ranks by target

    Args: 
        fpr (dataframe): a pandas dataframe with indices of probe names and a column corresponding to 
                        single probe ranks with column header 'rank'. may contain extra columns.
    Returns:
        target_fitness (dataframe): a pandas dataframe of target fitnesses weighted by rank.

    """

    df = fpr.copy()

    # get number of probes weights
    df.loc[:,'rank_weight'] = (df['rank'] - (n_probes_per_target + 1))**2

    # get probe weighted fitness
    df.loc[:,'weighted_fitness'] = df['fitness'] * df['rank_weight']

    # get target weighted fitness
    total_weight = sum([(i+1)**2 for i in range(n_probes_per_target)])
    target_fitness = df.groupby('target_id')['weighted_fitness'].sum()/total_weight

    return target_fitness

def weighted_target_pi( eij, fp, w0, n_probes_per_target=2, epsilon = 1e-6,
                        null_regex=".*NonTargeting.*",
                        library_file="~/dual_crispr/library_definitions/test_library_2_bpm.txt"):

    """Calculate the target level pi scores using a weighted rank based approach.
    TODO: this uses the ansatz of squared rank performance to calculate means, should we be doing this?

    Args:
        regex (str): a python regular expression to match null probes; default dual crispr style is '.*NonTargeting.*'
        fp (dataframe): a pandas dataframe with indices of probe names and a column corresponding to 
                        single probe fitnesses with column header 'fitness'. may contain extra columns.
        eij (dataframe): a pandas dataframe of size NxN with pi score estimates.
                        Columns and indices correspond to individual probes and the values of eij[i][j] is the 
                        dual probe interaction score. Matrix is symmetric about the diagonal.
        fc (dataframe): a pandas dataframe of size NxN with construct fitnesses.
                        Columns and indices correspond to individual probes and the values of fc[i][j] is the 
                        dual probe interaction score. Matrix is symmetric about the diagonal.
    Returns:
        eijm (dataframe): a pandas dataframe containing the pi scores
    """


    # get null mean
    null_mean = _get_null_dataframe(fp, null_regex = null_regex)['fitness'].mean()

    # substract the null mean from all the fitnesses
    fp['fitness'] = fp['fitness'] - null_mean

    # build probe to target lookup table
    library = _build_target_to_probe_table(library_file)

    # rank probes , 1=best, 2= second best .... 
    fpr = rank_probes(fp, library, null=True, null_regex = null_regex)

    # make boolean mask for expressed constructs
    expressed = w0>0

    # get construct weights
    construct_weights = ansatz_construct_weights(fpr, n_probes_per_target=n_probes_per_target)

    # calculated weighted pi scores, filtering on only those expressed
    eij_weighted = eij*construct_weights[expressed]

    # melt the weighted pi scores and add target ids
    eij_weighted_melt = _merge_in_target_ids(eij_weighted, library)

    # melt the construct_weights and merge in target ids
    construct_weights_melt = _merge_in_target_ids(construct_weights[expressed], library)

    # calculate sum of weighted pi scores by dual knock out target ids
    eij_sum = eij_weighted_melt.groupby(['target_id_1','target_id_2']).sum()

    # calculate normalizer from mean of target pairs weights 
    normalizer = construct_weights_melt.groupby(['target_id_1','target_id_2']).sum()

    # default the normalizer to some small value if ~0
    normalizer[normalizer['value'] < epsilon] = epsilon

    # normalize the constructs to the applied weights
    eijm = eij_sum / normalizer

    return eijm