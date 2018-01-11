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

def bootstrap(fc, pp, sdfc, w0, probes, 
    ag=2,
    tol=1e-3, 
    seed=None,
    maxiter=50,
    n_probes_per_target=2,
    epsilon = 1e-6,
    library_file="~/dual_crispr/library_definitions/test_library_2.txt",
    null_target_regex="NonTargetingControl",
    null_target_id="NonTargetingControl",
    null = True,
    verbose = False,
    testing = False,
    use_full_dataset_for_ranking = True
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

        fc_1 = np.array(r_fc(fc_0,sdfc_0, pp_0, seed))
        print(fc_1[:5,:5])
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

    if use_full_dataset_for_ranking:
        # get unweighted estimates for full dataset 
        # TODO: really shouldn't be running this on everybootstrap
        fp_0, fij_0, eij_0 = irls(fc, w0, ag=ag, tol=tol, maxiter=maxiter, verbose=verbose)
    else:
        fp_0 = None
        
    # get weighted pi scores and target fitness 
    pi_scores, target_fitness = weight_by_target(eij, fp, w0, probes,
                                                n_probes_per_target=n_probes_per_target,
                                                library_file = library_file,
                                                null_target_regex = null_target_regex,
                                                null_target_id = null_target_id,
                                                null = null,
                                                fp_0 = fp_0
                                                )
    # flatten for posterity
    pi_scores = pi_scores.stack().to_frame()

    return pi_scores, target_fitness

def run_iteration(fc, pp, sdfc, w0, probes,
    ag=2,
    tol=1e-3, 
    seed=None,
    maxiter=50,
    n_probes_per_target=2,
    epsilon = 1e-6,
    library_file="~/dual_crispr/library_definitions/test_library_2.txt",
    null_target_regex="NonTargetingControl",
    null_target_id="NonTargetingControl",
    null = True,
    niter=2,
    use_seed=False,
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
    pi_iter = None
    fitness_iter = None
    counter = 0
    while counter < niter:
        counter += 1

        if use_seed:
        	seed = counter
        else:
        	seed = None

        pi_scores, target_fitness = bootstrap(fc, pp, sdfc, w0, probes, 
                                                ag=ag, tol=tol, library_file=library_file, 
                                                seed=seed,
                                                maxiter=maxiter,
                                                n_probes_per_target=n_probes_per_target,
                                                epsilon = epsilon,
                                                null_target_regex=null_target_regex,
                                                null_target_id=null_target_id,
                                                null = null,
                                                verbose = verbose,
                                                testing= testing,
    						use_full_dataset_for_ranking = use_full_dataset_for_ranking
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


def summarize(pi_iter, fitness_iter, pp):

    from statsmodels.distributions.empirical_distribution import ECDF

    fitness_mean = fitness_iter.mean(axis=1)
    pi_mean = pi_iter.mean(axis=1)
    pi_iter_null = pi_iter.subtract(pi_mean, axis=0).values.flatten()

    enull = ECDF( np.concatenate((pi_iter_null, - pi_iter_null)) )
    emean = ECDF( pi_mean )

    fdr_left = np.minimum(1, enull(pi_mean)/emean(pi_mean))
    fdr_right = np.minimum(1, enull(-pi_mean)/(1-emean(pi_mean)))
    fdr = pd.DataFrame(list(zip(fdr_left, fdr_right)), index=pi_mean.index, columns=['fdr_left','fdr_right'])


    # get fdr for pairs
    fdr = compute_fdr(pi_iter)

    # calculate z score
    sd = pi_mean.sd()
    sd = sd.to_frame(columns=['std'])
    z = pi_mean.divide(sd)
    z = z.to_frame(columns=['Z'])

    # compute posterior probability
    pp = (pi_iter_null<pi_mean.abs()).mean(axis=1)
    
    # concatenate together
    x = pd.concat([pi_mean, fdr, sd, z, pp])

    # add genes as their own columns
    x.loc[:,'geneA'] = x.loc(axis=0)[,'target_id_1']
    x.loc[:,'geneB'] = x.loc(axis=0)[,'target_id_2']

    # merge in gene fitnesses
    x = pd.merge(x, fitness_mean, left_on='geneA', right_index=True, how="left")
    x = pd.merge(x, fitness_mean, left_on='geneB', right_index=True, how="left")

    return x


if __name__ == "__main__":

    # set up argument parser
    parser = argparse.ArgumentParser(usage = globals()['__doc__'])

    parser.add_argument("-v", "--verbose", action="store_true", default=False, \
                        help="Verbose output")
    parser.add_argument("--output", action="store", default=None, \
                        help="Directory to write results to.")
    parser.add_argument("--n_stds", action="store", default=2, type=float, \
                        help="Number of standard deviation to use in Tukey biweight normalization.")
    parser.add_argument("--tol", action="store", default=1e-3, type=float, \
                        help="Relative error tolerance")
    parser.add_argument("--max_irls_iter", action="store", default=50, type=int, \
                        help="Maximum IRLS iterations to perform")
    parser.add_argument("--n_probes_per_target", action="store", default=2, \
                        help="Maximum number of probes per target to use.") 
    parser.add_argument("--bootstrap_iterations", action="store", type=int, default=2, \
                        help="Number of bootstrap iterations to perform")
    parser.add_argument("--use_set_seed", action="store_true", default=False, \
                        help="Flag to use predefined seed for testing purposes.")
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




    # bootstrap weighted estimates
    pi_iter, fitness_iter = run_iteration(fc, pp, sdfc, w0, probes,
                                        ag=options.n_stds, tol=options.tol, 
                                        library_file=options.library_file,
                                        use_seed = options.use_set_seed,
    									n_probes_per_target = options.n_probes_per_target,
										niter=options.bootstrap_iterations,
										maxiter=options.max_irls_iter,
										verbose=options.verbose)

    mean_pi, mean_fitness = compute_mean(pi_iter, fitness_iter)
    if options.verbose:
        print("Pi scores:")
        print(mean_pi.head())
        print("Target Fitness:")
        print(mean_fitness.head())
