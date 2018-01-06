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

from irls import *
from weighted_pi import *

def bootstrap(fc, pp, sdfc, w0, probes, 
    ag=2, tol=1e-3, 
    seed=None, maxiter=50,
    n_probes_per_target=2,
    epsilon = 1e-6,
    library_file="~/dual_crispr/library_definitions/test_library_2.txt",
    null_target_regex="NonTargetingControl",
    null_target_id="NonTargetingControl",
    null = True,
    verbose = False
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

    # get unweighted estimates
    fp, fij, eij = irls(fc_1, w0, ag=ag, tol=tol, maxiter=maxiter, verbose=verbose)

    # get weighted pi scores and target fitness
    pi_scores, target_fitness = weight_by_target(eij, fp, w0, probes,
                                                n_probes_per_target=n_probes_per_target,
                                                library_file = library_file,
                                                null_target_regex = null_target_regex,
                                                null_target_id = null_target_id,
                                                null = null
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
    verbose=False
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
        	seed = None
        else:
        	seed = counter

        pi_scores, target_fitness = bootstrap(fc, pp, sdfc, w0, probes, 
                                                ag=ag, tol=tol, library_file=library_file, 
                                                seed=seed, maxiter=maxiter,
                                                n_probes_per_target=n_probes_per_target,
                                                epsilon = epsilon,
                                                null_target_regex=null_target_regex,
                                                null_target_id=null_target_id,
                                                null = null,
                                                verbose = verbose
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


def compute_mean(pi_iter, fitness_iter):
    """ Comput means across iterations
    """
    mean_pi = pi_iter.mean(axis=1)
    mean_fitness = fitness_iter.mean(axis=1)
    return mean_pi, mean_fitness


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