# """
# Docstring


# This is a rewrite of the iterative pi score bootstraping from Roman/Amanda's CTG pacakge.
# Ported to pure python.

# Todo:
#     method documentation
#     unit tests

# """
import re
import os
import argparse
import pandas as pd
import numpy as np

import rpy2.robjects.numpy2ri 
rpy2.robjects.numpy2ri.activate()

from irls import *
from weighted_pi import *

from statsmodels.distributions.empirical_distribution import ECDF

def _build_target_array(probes, names):
    """ Description
    """

    # get a list of all probes to targets from names datafame
    a = names[['probe_a_id','target_a_id']]
    b = names[['probe_b_id','target_b_id']]
    a.columns = ['probe_id','target_id']
    b.columns = ['probe_id','target_id']
    c = pd.concat([a,b],axis=0).drop_duplicates()

    # merge together and take target ids
    targets = np.array(pd.merge(pd.DataFrame(probes,columns=['probe_id']), c)['target_id'])
    return targets


def bootstrap(fc, pp, sdfc, w0, probes, targets,
    ag=2,
    tol=1e-3, 
    seed=None,
    maxiter=50,
    n_probes_per_target=2,
    epsilon = 1e-6,
    null_target_id="NonTargetingControl",
    null = True,
    verbose = False,
    testing = False,
    pre_computed_ranks = None
    ):
    """
    Perform 1 bootstrap of the pi score calculation based on posterior probabilities
    Args:
        fc (matrix):  A NxN matrix of observed fitnesses for each construct.
                        N is the number of probes.
                        Columns and rows are individual probes and the the value for [n,m] is the observed fitness of
                        the dual mutant comprised of probe m and probe n.  This fitness is the slope of a regression 
                        of log2 normalized abundances vs time.
                        Matrix is symmetric.
        w0 (matrix):  A NxN matrix of inital"weights" which is a boolean value to include a given construct pair in
                        the fit; 1 is True, include the pair, while 0 is False, do not include it.
                        Used to silence bad constructs.
                        N is the number of probes.
                        A probe pair or construct is given a 0 if it was not measured in the dataset or is deemed poor.
        probes (list):  A list of probe ids
        ag (int):       Which dimension to aggregate data across
        tol (float):    The error tolerance threshold which stops the iteration.
        maxit (int):    The maximum number of iterations to perform should 'tol' never be satisfied.
    """
    fc_0 = np.triu(fc)
    sdfc_0 = np.triu(sdfc)
    pp_0 = np.triu(pp)


    if testing:
        if seed is None:
            raise(BaseException("If testing for validation must provide a seed."))

        # basically use R random to mimick RS/A iterative pi scores
        # define r functions to get randoms with predefined seed
        rpy2.robjects.r('''
            gen_fc1 <- function(fc_0, sdfc_0, pp_0, iter) {
                utri<-upper.tri(fc_0)
                ntri<-sum(utri)
                nprobes = dim(fc_0)[1]
                fc_1<-matrix(0,nrow=nprobes,ncol=nprobes) 

                set.seed(iter)
                fc0<-fc_0[utri]+rnorm(ntri,sd=sdfc_0[utri]) 
                pp0<-pp_0[utri] 

                set.seed(iter)
                draw<-ifelse(runif(ntri)<pp0,1,0) 
                fc_1[utri]<-fc0*draw 
                fc_1<-fc_1+t(fc_1) 
                fc_1
            }
            ''')
        r_fc = rpy2.robjects.r['gen_fc1']
        fc_1 = np.array(r_fc(fc_0, sdfc_0, pp_0, seed))
    else:  
        # get random noise 
        if seed:
            np.random.seed(seed)
        noise = np.array([np.random.normal(loc=0, scale=sd) if sd>0 else 0 for sd in sdfc_0.flatten()]).reshape(sdfc_0.shape)

        # add noise to fc
        fc_0 = fc_0 + noise

        # decide whether to use base on posterior probability 
        # TODO: why randomly draw, should we set a confidence threshold
        if seed:
            np.random.seed(seed)

        include = pp_0 < np.random.rand(pp_0.shape[0], pp_0.shape[1])

        # multiply the construct matrix by boolean inclusion based on posterior proability
        fc_0 = fc_0 * include

        # make symmeteric
        fc_1 = fc_0 + fc_0.transpose()


    # get unweighted estimates using bootstrapped construct fitnesses
    fp, fij, eij = irls(fc_1, w0, ag=ag, tol=tol, maxiter=maxiter, verbose=verbose)


    # if use_full_dataset_for_ranking:
    #     # get unweighted estimates for full dataset 
    #     # TODO: really shouldn't be running this on everybootstrap
    #     fp_0, fij_0, eij_0 = irls(fc, w0, ag=ag, tol=tol, maxiter=maxiter, verbose=verbose)
    # else:
    #     fp_0 = None
        
    # get weighted pi scores and target fitness 
    pi_scores, target_fitness = weight_by_target(eij, fp, w0, probes, targets,
                                                n_probes_per_target=n_probes_per_target,
                                                null_target_id = null_target_id,
                                                null = null,
                                                pre_computed_ranks = pre_computed_ranks
                                                )

    # flatten for posterity
    pi_scores = pi_scores.stack().to_frame()

    return pi_scores, target_fitness

def run_iteration(fc, pp, sdfc, w0, probes, targets,
    ag=2,
    tol=1e-3, 
    maxiter=50,
    n_probes_per_target=2,
    epsilon = 1e-6,
    null_target_id="NonTargetingControl",
    null = True,
    niter=2,
    verbose=False,
    testing=False,
    use_full_dataset_for_ranking = True
    ):



    """
    Calculate the pi scores by iterative fitting

    Args: 
        fc (matrix): construct fitnesses
        pp (matrix): posterior probabilities 
        sdfc (matrix): the standard deviation of construct fitnesses 
        w0 (matrix): boolean matrix of good constructs
        niter (int): number of iterations to perform


    """

    if use_full_dataset_for_ranking:
        # if we want to use the full dataset for probe rankings, calculate it now
        # get initial fit of probe fitnesses based on construct fitnesses
        fp_0, fij_0, eij_0 = irls(fc, w0, ag=ag, tol=tol, maxiter=maxiter, verbose=verbose)

        # build index
        fp_index = pd.MultiIndex.from_tuples([(probes[i], targets[i]) for i in range(len(probes))],names=['probe_id','target_id']) 

        # make a dataframe out of the probe fitnesses
        fp_0 = pd.DataFrame(fp_0, index=fp_index, columns=['fitness'])
        
        if null:
            # get null mean
            null_mean = get_null_array(fp_0, null_target_id = null_target_id).mean(axis=0) 

            # substract the null mean from all the fitnesses
            fp_0 = fp_0 - null_mean

        # get ranks
        fpr_0 = rank_probes(fp_0, null=null, null_target_id= null_target_id)
    else:
        # if we don't want to use it then simply set the ranking datafame to None
        fpr_0 = None


    pi_iter = None
    fitness_iter = None
    counter = 0
    while counter < niter:
        counter += 1

        if verbose:
            print('Performing iteration: {}'.format(counter))
        if testing:
        	seed = counter
        else:
        	seed = None

        pi_scores, target_fitness = bootstrap(fc, pp, sdfc, w0, probes, targets,
                                                ag=ag,
                                                tol=tol,
                                                seed=seed,
                                                maxiter=maxiter,
                                                n_probes_per_target=n_probes_per_target,
                                                epsilon = epsilon,
                                                null_target_id=null_target_id,
                                                null = null,
                                                verbose = verbose,
                                                testing= testing,
                                                pre_computed_ranks = fpr_0
                                                )
        # add column labels corresponding to bootstrap
        pi_scores.columns = [counter]
        target_fitness.columns = [counter]

        # store results
        if pi_iter is None:
            pi_iter = pi_scores.copy()
            fitness_iter = target_fitness.copy()
        else:
            pi_iter = pd.concat([pi_iter, pi_scores], axis=1)
            fitness_iter = pd.concat([fitness_iter, target_fitness], axis=1)

    # return results of bootstrap
    return pi_iter, fitness_iter 

def compute_fdr(pi_iter):
    from statsmodels.distributions.empirical_distribution import ECDF

    pi_mean = pi_iter.mean(axis=1)
    pi_iter_null = pi_iter.subtract(pi_mean, axis=0).values.flatten()

    enull = ECDF( np.concatenate((pi_iter_null, - pi_iter_null)) )
    emean = ECDF( pi_mean )

    fdr_left = np.minimum(1, enull(pi_mean)/emean(pi_mean))
    fdr_right = np.minimum(1, enull(-pi_mean)/(1-emean(pi_mean)))
    fdr = pd.DataFrame(list(zip(fdr_left, fdr_right)), index=pi_mean.index, columns=['fdr_left','fdr_right'])

    return fdr

def summarize(pi_iter, fitness_iter):

    
    fitness_mean = fitness_iter.mean(axis=1).to_frame()
    pi_mean = pi_iter.mean(axis=1)
    pi_iter_null = pi_iter.subtract(pi_mean, axis=0)

    enull = ECDF( np.concatenate((pi_iter_null.values.flatten(), - pi_iter_null.values.flatten())) )
    emean = ECDF( pi_mean )

    fdr_left = np.minimum(1, enull(pi_mean)/emean(pi_mean))
    fdr_right = np.minimum(1, enull(-pi_mean)/(1-emean(pi_mean)))
    fdr = pd.DataFrame(list(zip(fdr_left, fdr_right)), index=pi_mean.index, columns=['fdr_left','fdr_right'])

    # get per pair sd
    sd = pi_iter.std(axis=1)

    # calculate z score
    pi_mean_sd = np.std(pi_mean)
    z = pi_mean.divide(pi_mean_sd)

    # compute posterior probability
    pp = pi_iter_null.abs().lt(pi_mean.abs(), axis=0).mean(axis=1)

    # concatenate together
    x = pd.concat([pi_mean, fdr, sd, z, pp], axis=1)
    x = x.reset_index()

    # merge in gene fitnesses
    x = pd.merge(x, fitness_mean, left_on='target_id_1', right_index=True, how="left")
    x = pd.merge(x, fitness_mean, left_on='target_id_2', right_index=True, how="left")
    x.columns = ['geneA','geneB','pi_mean', 'fdr_left', 'fdr_right', 'sd', 'z', 'pp', 'fA','fB']

    return x


