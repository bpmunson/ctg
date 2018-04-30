"""
Functions to calculated weighted mean of probe level fitness and pi scores
into target level estimates.
"""

import re
import os
import logging
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from collections import defaultdict
import scipy.sparse as sps



def rank_probes(fp, targets, null_target_id = "NonTargetingControl"):
    """ Rank probes by target according to fitness effects.

        Better probes are assumed to produce greater effect sizes 
        (either positive or negative) in fitness. However, for just the
        null target (scramble) the opposite is true: the best probes produce
        the smallest deviation from zero fitness.

        Args:
            fp (np.array): probe level fitness estimates
            targets (np.array): target names corresponding to probes
            null_target_id (str): name of the null/scramble target (gene)

        Returns:
            ranks (np.array): an array of integers corresponding to within 
                target based fitness ranks most likely 1.0, 2.0, 3.0, ... with 1 
                being the best 

    """
    # get unique targets
    uniq_targets = np.unique(targets)

    ranks = []
    # loop through each target, get the indicies, subset and rank
    for t in uniq_targets:
        ix = np.where(targets == t)[0]
        if t == null_target_id:
            # for null, best probes are closest to zero
            r = rankdata( abs( fp[ix]) ) 
        else:
            # for all other genes, best probes deviate the most from zero
            r = rankdata( - abs( fp[ix]) ) 
        # append ranks and index to growing list
        for i in range(len(r)):
            ranks.append((ix[i], r[i]))

    # convert to array
    ranks = np.array(ranks)
    # sort by index and take only ranks
    ranks = ranks[ranks[:,0].argsort()][:,1]

    return ranks
  
def mean_target_fitness(fp, ranks, targets, n_probes_per_target=2):
    """ Compute the weighted fitness per target, collapsing the probe based
        fitnesses according to the mean of the best probes
    
    Args: 
        fp (np.array): probe level fitness estimates
        ranks (np.array): an array of integers corresponding to within 
            target based fitness ranks most likely 1.0, 2.0, 3.0, ... with 1 
            being the best 
        targets (np.array): target names corresponding to probes
        n_probes_per_target (int): the desired number of probes to use per 
            target in calculating the weighted mean

    Returns:
        target_fitness (pd.dataframe): a pandas dataframe of target 
            fitnesses weighted by rank.

    """
    ranks_adj = n_probes_per_target + 1 - ranks
    # zero out negatives
    ranks_adj[ ranks_adj < 0 ] = 0
    # set all weights to 1
    rank_weight = np.ones(len(fp))
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
    return fitness

def mean_construct_weights(ranks, eij, n_probes_per_target=2):
    """Build the construct weight matrix given an array ranks by probe 

    Compute the mean weight for the construct array taking the best 
    n_probes_per_target for each target.  Here the mean is taken across all 
    accepted probes with no magnitude assigned to the rank beyond a binary
    inclusion.

    Args: 
        ranks (np.array): an array of integers corresponding to within 
            target based fitness ranks most likely 1.0, 2.0, 3.0, ... with 1 
            being the best
        eij (sp.sparse): probe level pi scores
        n_probes_per_target (int): the desired number of probes to use per 
                target in calculating the weighted mean
    Returns:
        construct_weights (ndarray): a NxN array of rank based weights for
            each probe pair.

    """
    # inverse ranks for weighting ... 
    # ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets
    # a weight of 4 or n_probes_per_target^2
    ranks_adj = n_probes_per_target + 1 - ranks
    # zero out negatives
    ranks_adj[ ranks_adj < 0 ] = 0
    # reset all the raks to 1
    ranks_adj[ ranks_adj>0 ] = 1
    nonzero = eij.nonzero()
    construct_weights = sps.csr_matrix(
        (
            ranks_adj[nonzero[0]] * ranks_adj[nonzero[1]],
            (nonzero[0], nonzero[1])
        ), shape=eij.shape)
    return construct_weights

def ansatz_target_fitness(fp, ranks, targets, n_probes_per_target=2):
    """ Compute the weighted fitness per target, collapsing the probe based
        fitnesses according to to their ranks
    
    Args: 
        fp (np.array): probe level fitness estimates
        ranks (np.array): an array of integers corresponding to within 
            target based fitness ranks most likely 1.0, 2.0, 3.0, ... with 1 
            being the best 
        targets (np.array): target names corresponding to probes
        n_probes_per_target (int): the desired number of probes to use per 
            target in calculating the weighted mean
    
    Returns:
        target_fitness (pd.dataframe): a pandas dataframe of target 
            fitnesses weighted by rank.

    """
    # inverse ranks for weighting ... 
    # ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets
    # a weight of 4 or n_probes_per_target^2
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
        # zero out the rank weights, if the probe based fitness is zero.
        for j in ix:
            if fpw[j] == 0:
                rank_weight[j] == 0
        target_fitness = fpw[ix].sum() / rank_weight[ix].sum()
        fitness.append(target_fitness)
    # convert to numpy array
    fitness = np.array(fitness)
    return fitness, uniq_targets

def ansatz_construct_weights(ranks, eij, n_probes_per_target=2):
    """Build the construct weight matrix given an array ranks by probe 

        This computes the outer cross product of the rank array adjusting for
        the number of probes we want to include for each target

    Args: 
        ranks (np.array): an array of integers corresponding to within 
            target based fitness ranks most likely 1.0, 2.0, 3.0, ... with 1 
            being the best
        eij (sp.sparse): probe level pi scores
        n_probes_per_target (int): the desired number of probes to use per 
                target in calculating the weighted mean

    Returns:
        construct_weights (ndarray): a NxN array of rank based weights for
            each probe pair.

    """
    # inverse ranks for weighting ... 
    # ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets
    # a weight of 4 or n_probes_per_target^2
    ranks_adj = n_probes_per_target + 1 - ranks
    # zero out negatives
    ranks_adj[ ranks_adj < 0 ] = 0
    nonzero = eij.nonzero()
    construct_weights = sps.csr_matrix(
        (
            ranks_adj[nonzero[0]] * ranks_adj[nonzero[1]],
            (nonzero[0], nonzero[1])
        ), shape=eij.shape)
    return construct_weights

def weight_by_target( eij, fp, w0, probes, targets,
    n_probes_per_target=2,
    epsilon = 1e-6,
    null_target_id="NonTargetingControl",
    pre_computed_ranks = None,
    method = "ansatz"
    ):
    """ Calculate weighted mean of probe fitnesses and pi scores at the 
        target/gene level.
    
        Args:
            eij (sp.sparse): probe level pi scores
            fp (np.array): probe level fitness estimates
            w0 (sp.sparse): initial construct weights matrix
            probes (np.array): probe names corresponding to fp/eij
            targets (np.array): target names corresponding to probes
            n_probes_per_target (int): the desired number of probes to use per 
                target in calculating the weighted mean
            epsilon (float): min value to add to weight to ensure no zero divide
            null_target_id (str): name of the null/scramble target (gene)
            pre_computed_ranks (np.array): An array of probe ranks to use 
                instead of computing them. These would likely imply use of the 
                ansatz method ranked probed for calculating the weighted mean.
                Ranks are expected to correspond to provided fp, eij, probes,
                and target arrays.
            method (str): the method to use in calculating weighted mean.
                default is ansatz or probe rankings.

        Returns:
            eijm (pd.DataFrame): pi scores by target-target pairs
            target_fitness (pd.DataFrame): target level fitness estimates

    """

    # subtract the null probes from all the fitnesses
    if null_target_id:
        # get the indicies of the null targets
        ix = np.where( targets == null_target_id)
        if len(ix)==0:
            raise AssertionError(
                'No null targets corresponding to "{}" were found'.format(
                    null_target_id))
        # get mean of the null probe fitnesses
        null_mean = fp[ix].mean()
        if np.isnan(null_mean):
            logging.error("Null mean not found. Ignoring.")
        else:
            # subtract off null mean from all fitnesses
            fp = fp - null_mean

    if pre_computed_ranks is None:
        # rank the probes 
        ranks = rank_probes(fp, targets, null_target_id = null_target_id)
    else:
        ranks = pre_computed_ranks

    if method=="mean":
        # compute the weighted target fitness
        target_fitness_values, target_fitness_labels = mean_target_fitness(
            fp, ranks, targets, n_probes_per_target=n_probes_per_target)

        # get construct weights
        construct_weights = mean_construct_weights(
            ranks, eij, n_probes_per_target=n_probes_per_target)
    else:
        # compute the weighted target fitness
        target_fitness_values, target_fitness_labels = ansatz_target_fitness(
            fp, ranks, targets, n_probes_per_target=n_probes_per_target)

        # get construct weights
        construct_weights = ansatz_construct_weights(
            ranks, eij, n_probes_per_target=n_probes_per_target)

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

    target_fitness = pd.DataFrame(target_fitness_values,
        index=target_fitness_labels)


    return eijm, target_fitness

