"""
Docstring


This is a rewrite of the pi score weighting from Roman/Amanda's CTG pacakge.
Ported to pure python.

Replicate fiteness correlation
dropping probes

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

from irls import *

from scipy.stats import rankdata
from collections import defaultdict

def magnitude_construct_weights(fpr, n_probes_per_construct=2):
    """ Use product of magnitude from zero as the weight """

    df = fpr.copy()

    # just like before have the option to zero out the worst performing probe for a given gene
    df.loc[:,'ranks_adj'] = n_probes_per_target + 1 - df.loc[:,'ranks']

    # zero out the negatives
    df['rank_adjusted'][ df['rank_adjusted'] < 0 ] = 0
    df['fitness'      ][ df['rank_adjusted'] < 0 ] = 0

    construct_weights = np.outer(df['fitness'], df['fitness'])

    return construct_weights

def magnitude_target_fitness(fpr, n_probes_per_target):
    """ Use the product of magnitude from zero as the weight """

    # make copy of input dataframe
    df = fpr.copy()

    # inverse ranks for weighting ... ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets a weight of 4 or n_probes_per_target^2
    df['rank_adjusted'] = (n_probes_per_target + 1 - df['rank'])

    # zero out negatives
    df['rank_adjusted'][ df['rank_adjusted'] < 0 ] = 0
    df['fitness'][ df['rank_adjusted'] < 0 ] = 0

    # square the ranks to get weights
    df['weighted_fitness'] = df['fitness']**2

    # multiply the probe fitness by the assigned rank weights
    df['weighted_fitness'] = df['fitness'] * df['rank_weight']

    # get target weighted fitness as the weighted mean of probe fitnesses
    target_fitness = df.groupby('target_id')['weighted_fitness'].sum() / df.groupby('target_id')['fitness'].sum()

    return target_fitness



def rank_probes(fp, targets, null=True, null_target_id = "NonTargetingControl"):
    """ Description

    Speed of this could be improved
    """
    # get unique targets
    uniq_targets = np.unique(targets)

    ranks = []
    # loop through each target, get the indicies, subset and rank
    for t in uniq_targets:
        ix = np.where(targets == t)[0]
        if t == null_target_id:
            r = rankdata( abs( fp[ix]) ) # for null, best probes are closest to zero
        else:
            r = rankdata( - abs( fp[ix]) ) # for all other genes, best probes deviate the most from zero
        # append ranks and index to growing list
        for i in range(len(r)):
            ranks.append((ix[i], r[i]))

    # convert to array
    ranks = np.array(ranks)
    # sort by index and take only ranks
    ranks = ranks[ranks[:,0].argsort()][:,1]

    return ranks
  

def ansatz_target_fitness(fp, ranks, targets, n_probes_per_target=2):
    """ Compute the weighted fitness per target, collapsing the probe based fitnesses according to to their ranks
    Args: 
        fpr (dataframe): an array of imputed probe base fitness values. with at least two columns 'fitness' and 'rank'
                     and a multiindex containing probe and target ids

    Returns:
        target_fitness (dataframe): a pandas dataframe of target fitnesses weighted by rank.

    """
    # inverse ranks for weighting ... ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets a weight of 4 or n_probes_per_target^2
    ranks_adj = n_probes_per_target + 1 - ranks
    # zero out negatives
    ranks_adj[ ranks_adj < 0 ] = 0
    # square the ranks to get weights
    rank_weight = ranks_adj**2
    # multiply the probe fitness by the assigned rank weights
    fpw = fp * rank_weight
    # get target weighted fitness as the weighted mean of probe fitnesses
    fitness = []
    uniq_targets = np.unique(targets)
    for t in np.unique(targets):
        ix = np.where(targets == t)[0]
        target_fitness = fpw[ix].sum() / rank_weight[ix].sum()
        fitness.append(target_fitness)
    # convert to numpy array
    fitness = np.array(fitness)
    return fitness, uniq_targets


def ansatz_construct_weights(ranks, eij, n_probes_per_target=2):
    """Build the construct weight matrix given an array ranks by probe 

    This computes the outer cross product of the rank array adjusting for the number of probes we want to include
    for each target

    Args: 
        ranks (np.array): an array of integers corresponding to within target based fitness ranks
                         most likely 1.0, 2.0, 3.0, ... with 1 being the best 
    Returns:
        construct_weights (ndarray): a NxN array of rank based weights for each probe pair.

    """
    # inverse ranks for weighting ... ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets a weight of 4 or n_probes_per_target^2
    ranks_adj = n_probes_per_target + 1 - ranks
    # zero out negatives
    ranks_adj[ ranks_adj < 0 ] = 0
    nonzero = eij.nonzero()
    construct_weights = sps.csr_matrix((ranks_adj[nonzero[0]] * ranks_adj[nonzero[1]], (nonzero[0], nonzero[1])), shape=eij.shape)
    return construct_weights


def weight_by_target( eij, fp, w0, probes, targets,
    n_probes_per_target=2,
    epsilon = 1e-6,
    null_target_id="NonTargetingControl",
    null = True,
    pre_computed_ranks = None
    ):
    """ Description
    """

    # subtract the null probes from all the fitnesses
    if null:
        # get the indicies of the null targets
        ix = np.where( targets == null_target_id)
        if len(ix)==0:
            raise AssertionError('No null targets corresponding to "{}" were found'.format(null_target_id))
        # get mean of the null probe fitnesses
        null_mean = fp[ix].mean()
        # subtract off null mean from all fitnesses
        fp = fp - null_mean

    if pre_computed_ranks is None:
        # rank the probes 
        ranks = rank_probes(fp, targets, null=null, null_target_id = null_target_id)
    else:
        ranks = pre_computed_ranks

    # compute the weighted target fitness
    target_fitness_values, target_fitness_labels = ansatz_target_fitness(fp, ranks, targets, n_probes_per_target=2)

    # get construct weights
    construct_weights = ansatz_construct_weights(ranks, eij, n_probes_per_target=n_probes_per_target)

    # make boolean mask for expressed constructs
    expressed = w0>0

    # calculated weighted pi scores, filtering on only those expressed
    construct_weights = construct_weights.multiply(expressed)
    eij = eij.multiply(construct_weights)


    # init storage container as nested defaultdict
    pi = defaultdict(float)
    weight = defaultdict(float)
    # convert to coordinate sparse matrix
    cx = sps.triu(eij).tocoo()
    # loop through nonzero elements
    for i,j,v in zip(cx.row, cx.col, cx.data):
        # get target ids
        target_a = targets[i]
        target_b = targets[j]
        # get construct weight 
        w = construct_weights[i,j]
        # store in dictionary
        pi[(target_a, target_b)]+=v
        weight[(target_a, target_b)]+=w

    # ensure weight normalizer is at least some minimum value
    for p in weight:
        if weight[p]<epsilon:
            weight[p] = epsilon

    # calculated weighted sum
    for p in pi:
        pi[p] = pi[p] / weight[p]

    # we currently expect a dataframe output, so make that for now
    eijm = pd.DataFrame(pi, index=['pi']).transpose()
    eijm.index.names = ['target_id_1', 'target_id_2']
    #eijm = eijm.reset_index().pivot(index='target_id_1', columns='target_id_2', values = 'pi')

    #target_fitness = pd.DataFrame(list(zip(target_fitness_labels, target_fitness_values)))
    target_fitness = pd.DataFrame(target_fitness_values, index=target_fitness_labels)


    return eijm, target_fitness

