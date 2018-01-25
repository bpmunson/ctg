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


def get_null_array( df, null_target_id = "NonTargetingControl"):
    """ Get the null array corresponding to rows of non targeting probes
    Args: 
        df (np.array): any array, n x m
        null_target_id (str): Name of null target, this was likely set to "NULL" in add_target_to_probe_index()
    Returns:
        null_df (np.array): returns a array QxM where Q is the number of null probes,
                            M is the original number of columns
    Raises:
        KeyError: if no rows corresponding to provided null_target_id were found
    """
    if (df.index.get_level_values('target_id')==null_target_id).any():
        null_df = df.xs(null_target_id, level="target_id")
    else:
        raise(KeyError("No null targets found corresponding to '{}'".format(null_target_id)))
    return null_df

def rank_probes(fp, null=True, null_target_id = "NonTargetingControl"):
    """Rank the probes by deviation from zero
    Todo:
    should we require the library_definition?
    add special case for null probes where we want the closest ones to zero
    should the null_target be a list of nulls?

    Args:
        fp (dataframe): an array of imputed probe base fitness values. with at least one columns 'fitness' 
                     and a multiindex containing probe and target ids
        null (bool): treat the null probes differently
        null_target_id (str): a string corresponding to which sample was the null sample

    Returns:
        fpr (dataframe): the original fp dataframe with an additional column corresponding to probe rnaks per target
    """

    # make a copy of the input dataframe
    fpr = fp.copy()

    # get absolute value of fitness
    fpr.loc[:,'fitness_abs'] = fpr.loc[:,'fitness'].abs()

    # get rank
    # for the knockout probesthe assumption here is that we want the probes which deviate furthest from zero 
    fpr.loc[:,'rank'] = fpr.groupby('target_id')['fitness_abs'].rank(method='max', ascending=False)

    # do something special for the null probes
    if null:
        # replace the ranks for the null probes
        # for the null probes we want to take the probes which are closest to zero 
        idx = pd.IndexSlice
        
        # get only null target rows
        if (fpr.index.get_level_values('target_id')==null_target_id).any():
            tmp_a = fpr.loc(axis=0)[idx[:, null_target_id]]
        else:
            raise(KeyError("No null targets found corresponding to '{}'".format(null_target_id)))

        # reverse the rank sorting for null probes, ie closer to zero is better
        tmp_b = tmp_a.groupby('target_id')['fitness_abs'].rank(method='min', ascending=True)
        # replace values
        fpr.loc[idx[:, null_target_id],'rank'] = tmp_b 

    # remove absolute column 
    fpr = fpr.drop(labels ='fitness_abs', axis=1)

    return fpr

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

    # inverse ranks for weighting ... ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets a weight of 4 or n_probes_per_target^2
    ranks_adj = n_probes_per_target + 1 - ranks

    # zero out negatives
    ranks_adj[ ranks_adj < 0 ] = 0

    # calculate outer cross product to get construct weights
    construct_weights = np.outer(ranks_adj, ranks_adj)

    return construct_weights

def ansatz_target_fitness(fpr, n_probes_per_target=2):
    """ Compute the weighted fitness per target, collapsing the probe based fitnesses according to to their ranks
    Args: 
        fpr (dataframe): an array of imputed probe base fitness values. with at least two columns 'fitness' and 'rank'
                     and a multiindex containing probe and target ids

    Returns:
        target_fitness (dataframe): a pandas dataframe of target fitnesses weighted by rank.

    """

    # make copy of input dataframe
    df = fpr.copy()

    # inverse ranks for weighting ... ie best probe gets a rank of n_probes_per_target, normally 2 
    # so the best pair of probes for each target pair gets a weight of 4 or n_probes_per_target^2
    df['rank_adjusted'] = (n_probes_per_target + 1 - df['rank'])

    # zero out negatives
    df['rank_adjusted'][ df['rank_adjusted'] < 0 ] = 0

    # square the ranks to get weights
    df['rank_weight'] = df['rank_adjusted']**2

    # multiply the probe fitness by the assigned rank weights
    df['weighted_fitness'] = df['fitness'] * df['rank_weight']

    # get target weighted fitness as the weighted mean of probe fitnesses
    target_fitness = df.groupby('target_id')['weighted_fitness'].sum() / df.groupby('target_id')['rank_weight'].sum()

    return target_fitness

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

def weight_by_target( eij, fp, w0, probes, targets,
    n_probes_per_target=2,
    epsilon = 1e-6,
    null_target_id="NonTargetingControl",
    null = True,
    fp_0 = None
    ):

    """ Calculated the weighted average of pi scores and fitness 
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
        eijm (pandas dataframe): a dataframe containing the pi scores in square format
        target_fitnesses (pandas series): a named series containing gene/target fitnesses 
    """

    # make index out of targets and probes
    fp_index = pd.MultiIndex.from_tuples([(probes[i], targets[i]) for i in range(len(probes))],names=['probe_id','target_id'])

    # make a dataframe out of the probe fitnesses
    fp = pd.DataFrame(fp, index=fp_index, columns=['fitness'])

    if null:
        # get null mean
        null_mean = get_null_array(fp, null_target_id = null_target_id).mean(axis=0) 

        # substract the null mean from all the fitnesses
        fp = fp - null_mean

    # use the precompute ranks instead of those calculated here
    if fp_0 is not None:
        # TODO: clean up this repeat of code, don't recalculate the rankings everytime?
        fp_0 = pd.DataFrame(fp_0, index=fp_index, columns=['fitness'])

        if null:
            # get null mean
            null_mean = get_null_array(fp_0, null_target_id = null_target_id).mean(axis=0) 

            # substract the null mean from all the fitnesses
            fp_0 = fp_0 - null_mean

        # rank probes , 1=best, 2= second best .... 
        fpr = rank_probes(fp_0, null=null, null_target_id= null_target_id)

        # reassign fitness column names
        fpr.columns = ['fitness_0', 'rank']

        # concatenate and drop the full fitness column
        fpr = pd.concat([fp, fpr], axis=1)
        fpr = fpr.drop('fitness_0', axis=1)
    else:
        # rank probes , 1=best, 2= second best .... 
        fpr = rank_probes(fp, null=null, null_target_id= null_target_id)

    # compute the weighted target fitness
    target_fitness = ansatz_target_fitness(fpr, n_probes_per_target=2)

    # get construct weights
    construct_weights = ansatz_construct_weights(fpr['rank'], n_probes_per_target=n_probes_per_target)

    # make copies of the index and rename levels
    rows = fp.index.copy()
    rows.names = ['probe_id_1','target_id_1']
    columns = fp.index.copy()
    columns.names = ['probe_id_2','target_id_2']

    # make a dataframe out of the construct weights
    construct_weights = pd.DataFrame(construct_weights, columns=columns, index=rows )

    # make a dataframe out of the pi scores
    eij = pd.DataFrame(eij,  columns=columns, index=rows )
            
    # make boolean mask for expressed constructs
    expressed = w0>0
    expressed = pd.DataFrame(expressed,  columns=columns, index=rows )

    # calculated weighted pi scores, filtering on only those expressed
    eij_weighted = eij*construct_weights[expressed]

    # calculate sum of weighted pi scores by dual knock out target ids
    eij_sum = eij_weighted.stack(level=eij_weighted.columns.names)\
                    .groupby(['target_id_1','target_id_2']).sum()

    # calculate normalizer from mean of target pairs weights 
    normalizer = construct_weights[expressed].stack(level=construct_weights.columns.names)\
                    .groupby(['target_id_1','target_id_2']).sum()

    # default the normalizer to some small value if ~0 to zero out any non weigthed pi scores
    normalizer[normalizer < epsilon] = epsilon

    # normalize the constructs to the applied weights
    eijm = eij_sum / normalizer

    # cast pi scores into matrix (wide) format 
    eijm = eijm.to_frame()
    eijm.columns = ["pi"]
    eijm = eijm.reset_index().pivot(index='target_id_1', columns='target_id_2', values = 'pi')

    # finally return results
    return eijm, target_fitness



if __name__ == "__main__":

    # set up argument parser
    parser = argparse.ArgumentParser(usage = globals()['__doc__'])

    parser.add_argument("-v", "--verbose", action="store_true", default=False, \
                        help="Verbose output")
    parser.add_argument("-f", "--construct_fitness", action="store", default=None, required=True, \
                        help="Path to construct fitness matrix.")
    parser.add_argument("-s", "--construct_fitness_std", action="store", default=None, required=True, \
                        help="Path to construct fitness standard deviation.")
    parser.add_argument("-p", "--posterior_probability", action="store", default=None, required=True, \
                        help="Path to pi score posteriors.")
    parser.add_argument("-w", "--construct_weights", action="store", default=None, required=True, \
                        help="Path to construct weights.")
    parser.add_argument("-l", "--library_file", action="store", required=True, \
                        default="~/dual_crispr/library_definitions/test_library_2.txt", \
                        help="Path to library defintion.")
    parser.add_argument("-o", "--output", action="store", default=None, \
                        help="Directory to write results to.")
    parser.add_argument("-a", "--n_stds", action="store", default=2, \
                        help="Number of standard deviation to use in Tukey biweight normalization.")
    parser.add_argument("--tol", action="store", default=1e-3, \
                        help="Relative error tolerance")
    parser.add_argument("--maxiter", action="store", default=50, \
                        help="Maximum IRLS iterations to perform")
    parser.add_argument("-n", "--n_probes_per_target", action="store", default=2, \
                        help="Maximum number of probes per target to use.")

    # parse arguments
    options = parser.parse_args()

    if not options.construct_fitness and not options.construct_weights:
        raise BaseException("No input files provided.")

    # load input files
    fc = load_matrix_csv(options.construct_fitness)
    w0 = load_matrix_csv(options.construct_weights)
    sdfc = load_matrix_csv(options.construct_fitness_std)
    pp = load_matrix_csv(options.posterior_probability)

    # TODO: functionalize the probe label loadings
    probes = list(pd.read_csv(options.construct_fitness, header=0, index_col=0).columns)


    # compute initial fitness and interaction estimates
    fp, fij, eij = irls(fc, w0,
                        ag=options.n_stds,
                        tol=options.tol,
                        maxiter=options.maxiter,
                        verbose=options.verbose)


    # compute weighted estimates
    pi_scores, target_fitness = weight_by_target(eij, fp, w0, probes,
                                                n_probes_per_target=options.n_probes_per_target,
                                                library_file = options.library_file,
                                                null_target_regex = "NonTargeting",
                                                null_target_id = "0",
                                                null = True
                                                )



    if options.output:
        pi_scores.to_csv(os.path.join(options.output, "TestSet8_pi_scores.csv"), sep=",")
        target_fitness.to_csv(os.path.join(options.output, "TestSet8_f_target.csv"), sep=",")



