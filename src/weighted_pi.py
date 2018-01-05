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


def load_library_file(fp, sep="\t", index_col=0, header=0, comment="#"):
    """ Loads libary definition from a csv file

        TODO: assert this is a good library file, ie at a very mimumum validate the expected header

    """
    library = pd.read_csv(fp, sep=sep, index_col=index_col, header = header, comment=comment)
    return library

def get_null_array( df, null_target_id = "NULL"):
    """ Get the null array corresponding to rows of non targeting probes
    Args: 
        df (np.array): any array, n x m
        null_target_id (str): Name of null target, this was likely set to "NULL" in add_target_to_probe_index()
    Returns:
        null_df (np.array): returns a array QxM where Q is the number of null probes,
                            M is the original number of columns
    Raises:
        BaseException: if not rows corresponding to provided null_target_id were found
    """
    if (df.index.get_level_values('target_id')==null_target_id).any():
        null_df = df.xs(null_target_id, level="target_id")
    else:
        raise(KeyError("No null targets found corresponding to '{}'".format(null_target_id)))
    return null_df

def add_target_to_probe_index(  index,  
                                library=None,
                                library_file=None, \
                                null=True,
                                null_target_regex="NonTargetingControl",
                                null_target_id="NULL" ):
    """Description
    """
    if library is None:
        if library_file is None:
            raise BaseException("Must supply either a libary dataframe or a path to the library definition csv.")

        library = load_library_file(library_file)

    # TODO: remove this dependenacy of a 0 start from the input files
    #index = pd.Index([re.sub("0NonTargeting", "NonTargeting", i) for i in index] )

    # make a dataframe from the index
    i = pd.DataFrame(index, columns=['idx'])

    # get only relvant columns from the library definition
    libr_a = library[['probe_a_id','target_a_id']]
    libr_a.columns = ['probe_id','target_id']
    libr_b = library[['probe_b_id','target_b_id']]
    libr_b.columns = ['probe_id','target_id']
    libr = pd.concat([libr_a, libr_b])
    libr = libr.reset_index(drop=True)
    lib = libr.drop_duplicates()

    if null:
        # try and replace null targets
        # TODO: this raises a copy slice warning...
        lib['target_id'].replace(".*{}.*".format(null_target_regex), null_target_id, regex=True, inplace=True)
        # check to see if we found any
        if not (lib['target_id'] == null_target_id ).any():
            raise(KeyError("No null probes found corresponding to '{}'.".format(null_target_regex)))

    # merge the two together on probe_id
    res = pd.merge(i, lib, left_on="idx", right_on="probe_id", how="left")

    # set the index on the merge datfrmae 
    res.set_index(['probe_id','target_id'], inplace = True)

    # return the new index
    return res.index

def rank_probes(fp, null=True, null_target_id = "NULL"):
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
            raise(BaseException("No null targets found corresponding to '{}'".format(null_target_id)))

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

def weight_by_target( eij, fp, w0, probes,
    n_probes_per_target=2,
    epsilon = 1e-6,
    library_file="~/dual_crispr/library_definitions/test_library_2.txt",
    null_target_regex="NonTargetingControl",
    null_target_id="NULL",
    null = True):

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
        eijm (np.array): a matrix containing the pi scores
        target_fitnesses 
    """

    # make a dataframe out of the probe fitnesses
    fp = pd.DataFrame(fp, index=probes, columns=['fitness'])

    # get the lbrary definitions file 
    library = load_library_file(library_file, sep="\t", header=0, index_col=0, comment="#")

    # add target id to probe fitness dataframe
    # while change the NonTargeting names to NULL
    fp.index = add_target_to_probe_index(   fp.index,
                                            library=library,
                                            null=null,
                                            null_target_id = null_target_id,
                                            null_target_regex=null_target_regex)

    # make copies of the index and rename levels
    rows = fp.index.copy()
    rows.names = ['probe_id_1','target_id_1']
    columns = fp.index.copy()
    columns.names = ['probe_id_2','target_id_2']

    if null:
        # get null mean
        null_mean = get_null_array(fp, null_target_id = null_target_id).mean(axis=0) 

        # substract the null mean from all the fitnesses
        fp = fp - null_mean

    # rank probes , 1=best, 2= second best .... 
    fpr = rank_probes(fp, null=null, null_target_id= null_target_id)

    # compute the weighted target fitness
    target_fitness = ansatz_target_fitness(fpr, n_probes_per_target=2)

    # get construct weights
    construct_weights = ansatz_construct_weights(fpr['rank'], n_probes_per_target=n_probes_per_target)

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
    normalizer = construct_weights.stack(level=construct_weights.columns.names)\
                    .groupby(['target_id_1','target_id_2']).sum()

    # default the normalizer to some small value if ~0 to zero out any non weigthed pi scores
    normalizer[normalizer < epsilon] = epsilon

    # normalize the constructs to the applied weights
    eijm = eij_sum / normalizer

    # finally return results
    return eijm, target_fitness



if __name__ == "__main__":

    # set up argument parser
    parser = argparse.ArgumentParser(usage = globals()['__doc__'])
    parser.add_argument("-v", "--verbose", action="store_true", default=False, \
                        help="Verbose output")
    parser.add_argument("-f", "--construct_fitness", action="store", default=None, \
                        help="Path to construct fitness.")
    parser.add_argument("-s", "--construct_fitness_std", action="store", default=None, \
                        help="Path to construct fitness standard deviation.")
    parser.add_argument("-p", "--posterior_probability", action="store", default=None, \
                        help="Path to pi score posteriors.")
    parser.add_argument("-w", "--construct_weights", action="store", default=None, \
                        help="Path to construct weights.")
    parser.add_argument("-l", "--library_file", action="store", \
                        default="~/dual_crispr/library_definitions/test_library_2.txt", \
                        help="Path to library defintion.")
    parser.add_argument("-o", "--output", action="store", default=None, \
                        help="Directory to write results to.")
    parser.add_argument("-a", "--n_stds", action="store", default=2, \
                        help="Number of standard deviation to use in Tukey biweight normalization")
    parser.add_argument("--tol", action="store", default=1e-3, \
                        help="Relative error tolerance")
    parser.add_argument("--maxiter", action="store", default=50, \
                        help="Maximum IRLS iterations to perform")
    parser.add_argument("--bootstrap_iterations", action="store", default=2, \
                        help="Number of bootstrap iterations to perform")
    # parse arguments
    options = parser.parse_args()

    if not options.construct_fitness and not options.construct_weights:
        raise BaseException("No input files provided.")

    # load input files
    fc = load_construct_fitnesses(options.construct_fitness)
    w0 = load_initial_weights(options.construct_weights)
    sdfc = load_construct_fitnesses(options.construct_fitness_std)
    pp = load_construct_fitnesses(options.posterior_probability)

    # TODO: functionalize the probe label loadings
    probes = list(pd.read_csv(options.construct_fitness, header=0, index_col=0).columns)

    library_file = options.library_file
    epsilon=1e-6
    n_probes_per_target=2

    # compute initial fitness and interaction estimates
    fp, fij, eij = irls(fc, w0,
                        ag=options.n_stds,
                        tol=options.tol,
                        maxiter=options.maxiter,
                        verbose=options.verbose)


    # compute weighted estimates
    pi_scores, target_fitness = weight_by_target(eij, fp, w0, probes, library_file = options.library_file,
                                                null_target_regex = "NonTargeting",
                                                null_target_id = "NULL"
                                                )



