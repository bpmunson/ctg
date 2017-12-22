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



def _get_null_array(df, probes, indices = None, null_regex=".*NonTargeting.*"):
    """ Get the null array corresponding to rows of non targeting probes
    Args: 
        df (np.array): any array, n x m
        indices (list): True to return indices of null probes in df,
                        False to return probe names of null probes in df
    Returns:
        null_df (np.array): returns a array OxM where O is the number of null probes,
                            M is the original number of columns
    """
    
    if indices is None:
        # if no index range was included use regex 
        indices = _get_null_indices(probes, null_regex = null_regex)

    null_df = df[indices]

    return null_df

def _get_null_indices(probes, null_regex=".*NonTargeting.*"):
    """ Get the null indices of the dataframe corresponding to non targeting probes
    Args:
        probes (list): the probe names for the entire library
        regex (str): a python regular expression to match null probes; default dual crispr style is '.*NonTargeting.*'
    Returns:
        ids (list): a list of indices corresponding to the null target in the probe list
    """

    r = re.compile(null_regex)

    ids = [i for i, probe in enumerate(probes) if r.match(probe) ]

    return ids

def _get_null_probes_names(probes, null_regex=".*NonTargeting.*"):
    """ Get the null probe names of the full probe list corresponding to non targeting probes
    Args:
        probes (list): the probe names for the entire library
        regex (str): a python regular expression to match null probes; default dual crispr style is '.*NonTargeting.*'
    Returns:
        names (list): a list of probe names corresponding to the null target
    """

    r = re.compile(null_regex)
    names = [probe for probe in probes if r.match(probe)]
    
    return names

def _get_target_to_probe_pairs(fn, sep="\t", header=0, comment = "#", null_regex=".*NonTargeting.*"):
    """ Loads library definition
        Todo: this assumes a CTG style library definition, should probably abstract it
    Args:
        fn (str): path to library definition file
    Returns:
        probe_to_target (list): a list of (probe, target) tuples
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
    #c.drop('probe_id', axis=1, inplace=True)

    # replace NonTargetingNNNN target id 
    null_probes = _get_null_probes_names(c['probe_id'], null_regex = null_regex )
    c.loc[null_probes, 'target_id'] = "0NonTargeting" # TODO: change this hardcoding

    probe_to_target = c['probe_id','target_id'].values

    return probe_to_target

def rank_probes(fp, probes, probe_to_target, null=True, null_probes = None, null_regex = ".*NonTargeting.*"):
    """Rank the probes by 
    Todo:
    should we require the library_definition?
    add special case for null probes where we want the closest ones to zero
    Args:
        fp (np.array): an array of imputed probe base fitness values.
        probes (list): the probe names for the entire library
        probe_to_target (list): a list of (probe, target) tuples
        null (bool): treat the null probes differently
        null_probes (list): a list of probe names corresponding to the null target
        null_regex (str): a regex search string to find the null probe names in probe list 

    Returns:
        ranks (np.array): an array of integers corresponding to within target based fitness ranks
                         most likely 1.0, 2.0, 3.0, ... with 1 being the best 
    """
    
    # make probe fitnesses into dataframe with index being probe ids
    fp_df = pd.DataFrame(fp, index=probes, columns='fitness')

    # make target to probe pairs a dataframe with index being probe ids
    library = pd.DataFrame(probe_to_target, columns = ["probe_id","target_id"])
    library.index = library['probe_id']

    # add target ids to the fitness dataframe
    fpc = pd.merge(fp_df, library, left_index = True, right_index=True)

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
 
    ranks = np.array(fpc['rank'])

    return ranks

def ansatz_construct_weights(ranks, n_probes_per_target=2):
    """Build the construct weight matrix given an array ranks by probe 

    This computes the outer cross product of the rank array adjusting for the number of probes we want to include
    for each target

    Args: 
        ranks (np.array): an array of integers corresponding to within target based fitness ranks
                         most likely 1.0, 2.0, 3.0, ... with 1 being the best 
    Returns:
        construct_weights (ndarray): a NxN array of rank based weights for each probe pair.

    """

    ranks_adj = ranks - (n_probes_per_target + 1)

    construct_weights = np.outer(ranks_adj, ranks_adj)

    # construct_weights = np.zeros((n,n))
    # #rank_iloc = list(fpr.columns).index('rank')
    # for i in range(n):
    #     for j in range(n):
    #         # so
    #         # rank(1) * rank(1) is the best and gets a weight of 4
    #         # rank(2) * rank(1) is second best and gets a weight of 2
    #         # rank(2) * rank(2) is next and gets a weight of 1
    #         # rank(3) * rank(3) is worst and we wont use it so it gets a rank of 0

    #         # abstract this for any number of probes per construct
    #         construct_weights[i][j] =   (ranks[i] - (n_probes_per_target+1)) * \
    #                                     (ranks[j] - (n_probes_per_target+1))
    return construct_weights

def ansatz_target_fitness(fp, probes, ranks, n_probes_per_target=2):
    """ Compute the weighted fitness per target, collapsing the probe based fitnesses according to to their ranks

    Args: 
        fp (np.array): an array of imputed probe base fitness values.
        ranks (np.array): an array of integers corresponding to within target based fitness ranks
                         most likely 1.0, 2.0, 3.0, ... with 1 being the best 

    Returns:
        target_fitness (dataframe): a pandas dataframe of target fitnesses weighted by rank.

    """

    # make probe fitnesses into dataframe with index being probe ids
    df = pd.DataFrame(zip(fp, ranks), index=probes, columns=['fitness','probes'])

    # make target to probe pairs a dataframe with index being probe ids
    library = pd.DataFrame(probe_to_target, columns = ["probe_id","target_id"])
    library.index = library['probe_id']

    # add target ids to the fitness dataframe
    df = pd.merge(df, library, left_index = True, right_index=True)

    # get number of probes weights
    df.loc[:,'rank_weight'] = (df['rank'] - (n_probes_per_target + 1))**2

    # get probe weighted fitness
    df.loc[:,'weighted_fitness'] = df['fitness'] * df['rank_weight']

    # get target weighted fitnes
    target_fitness = df.groupby('target_id')['weighted_fitness'].sum()/df.groupby('target_id')['rank_weight'].sum()

    return target_fitness

def weighted_target_pi( eij, fp, w0, probes, n_probes_per_target=2, epsilon = 1e-6,
    null_regex=".*NonTargeting.*",
    library_file="~/dual_crispr/library_definitions/test_library_2_bpm.txt"):

    """Calculate the target level pi scores using a weighted rank based approach.
    TODO: this uses the ansatz of squared rank performance to calculate means, should we be doing this?

    Args:
        regex (str): a python regular expression to match null probes; default dual crispr style is '.*NonTargeting.*'
        fp (np.array): an array of imputed probe base fitness values.
        eij (np.array): a pandas dataframe of size NxN with pi score estimates.
                        Columns and indices correspond to individual probes and the values of eij[i][j] is the 
                        dual probe interaction score. Matrix is symmetric about the diagonal.
        fc (np.array): a pandas dataframe of size NxN with construct fitnesses.
                        Columns and indices correspond to individual probes and the values of fc[i][j] is the 
                        dual probe interaction score. Matrix is symmetric about the diagonal.
        probes (list): the probe names for the corresponding arrays
    Returns:
        eijm (np.array): a matrix containing the pi scores
    """


    # get null mean
    null_mean = _get_null_array(fp, null_regex = null_regex).mean()

    # substract the null mean from all the fitnesses
    fp = fp - null_mean

    # build probe to target lookup table
    probe_to_target = _get_target_to_probe_pairs(library_file)

    # rank probes , 1=best, 2= second best .... 
    ranks = rank_probes(fp, probes, probe_to_target, null=True, null_regex = null_regex)

    # make boolean mask for expressed constructs
    expressed = w0>0

    # get construct weights
    construct_weights = ansatz_construct_weights(ranks, n_probes_per_target=n_probes_per_target)

    # make a dataframe out of the pi scores
    eij = pd.DataFrame(eij, columns = probes, index=probes)
    
    # calculated weighted pi scores, filtering on only those expressed
    eij_weighted = np.ma.masked_array(eij*construct_weights, mask=expressed)

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

def solve(niter):
    """
    Calculate the pi scores by iterative fitting

    Args: 
        fc (matrix): construct fitnesses
        pp (matrix): posterior probabilities 
        sdfc (matrix): the standard deviation of construct fitnesses 
        w0 (matrix): boolean matrix of good constructs
        niter (int): number of iterations to perform


    """

    # get only upper triangles 
    fc_0 = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
    sdfc_0 = sdfc.where(np.triu(np.ones(sdfc.shape)).astype(np.bool))
    pp_0 = pp.where(np.triu(np.ones(sdfc.shape)).astype(np.bool))

    # get random noise 
    noise = np.array([np.random.normal(loc=0, scale=sd) if sd>0 else 0 for sd in sdfc_0.values.flatten()]).reshape(sdfc_0.shape)

    # add noise to fc
    fc_0 = fc_0 + noise

    # decide whether to use base on posterior probability 
    # TODO: why randomly draw, should we set a confidence threshold
    include = pp_0 < np.random.rand(pp_0.shape[0], pp_0.shape[1])

    # multiply the construct matrix by boolean inclusion based on posterior proability
    fc_0 = fc_0 * include

    # make symmeteric
    fc_1 = fc_0 + tc.transpose()

    # get new estimates
    fp, fij, eij = irls(fc_1, w0)

    # get weighted target level fitnesses
    target_fitnesses = ansatz_target_fitness(fp, library, n_probes_per_target=2)

    # get weighted pi scores
    pi_scores = weighted_target_pi(eij, fp, w0)
    
    return target_fitness, pi_scores



