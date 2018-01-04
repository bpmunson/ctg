"""
Docstring


This is a rewrite of the irls function from Roman/Amanda's CTG pacakge.
Ported to pure python.

Todo:
    method documentation
    verify biweight
    unit tests
    change ag to float but keep options for number of standard devaitions away
"""

import os
import argparse
import pandas as pd
import numpy as np


def load_construct_fitnesses(fp, sep=",", index_col=0, header=0):
    """ Loads precomputed constuct fitnesses from a csv file
    """
    
    fc = np.loadtxt(fp, delimiter=sep)
    return fc

def load_initial_weights(fp, sep=",", index_col=0, header=0):
    """ Loads precomputed weights for construct inclusion in fitting.
    """
    
    w0 = np.loadtxt(fp, delimiter=sep)
    return w0

def filter_weights(w, w0):
    """Filter weights for bad constructs and constructs with the only one target
    """
    
    # replace the diagonal (ie same gene on both probes ) with zero ... or do not include
    np.fill_diagonal(w, 0)

    # zero out previous determined bad constructs 
    w = w * w0

    return w

def biweight(eij, ag=2, expressed=None):
    """ Normalize weights according to a modified Tukey's Biweight
    """
    if expressed is None:
        # create a mask of all true
        expressed = np.ones(eij.shape).astype(np.bool)
        
    # calculate standard deviation of pi scores (residuals) for expressed constructs
    sd = eij[expressed].std()
    
    # calculate new weights, aparently using a modified biweight model 
    # TODO, why not multiply by eij in front?  
    w = (1-eij**2/(ag*sd)**2)**2
    
    # zero out anythin more than ag standard deviations from zero
    w[abs(eij)>(ag*sd)]=0

    return w



def construct_system_ab(fc, w, err=1e-6):
    """ Construct the system of equations to solve
    """
    n = fc.shape[0]
    # initialize the imputed weights matrix
    A = w.copy()

    # replace diagnol with total number of constructs used for the given probe,
    # plus some small error term
    for i in range(n):
        A[i,i] = w[:,i].sum()+err

    # get new construct fitnesses to fit against
    # initialize a vector with sum of all construct fitnesses for every probe
    # if the assumption of a normally distributed genetic interactions about zero is met 
    # this will be 0 or close to it ... TODO: should it be the mean? 
    b = np.array((fc*w).sum(axis=0))
    
    return A, b

def solve_iteration(A, b):
    """ 
    Solves for the individual probe fitnesses
    """

    # find single probe fitnesses which satisfy expectation, you additive property
    x = np.linalg.solve(A, b)

    # init probe estimates
    n = A.shape[0]
    fij = np.zeros((n,n))

    # fill expected fitness and pi score matrices
    for i in range(n):
        for j in range(i+1,n):
            fij[i,j] = x[i]+x[j]

    # convert to a pandas data frame
    #fij = pd.DataFrame(fij, columns = A.columns, index = A.index)

        
    # make the expected fitness matrix symmetric
    fij = fij + np.transpose(fij)
    eij = fc - fij

    return x, fij, eij

def irls(fc, w0, ag=2, probes=None, tol=1e-3, maxiter=50, verbose=False):
    """ The iterative least squares fitting function of single gene fitnesses

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

    Returns:
        fp (list):      The individual probe fitnesses. A list of floats of length N, the number of probes.
        pi (matrix)     A symmetric matrix of floats corresponding to the raw pi scores resulting from the fit.  
                        Size if NxN where N is the number of probes.
        fij (matrix)    A symmetric matrix of floats corresponding to the expected fitnesses for each contruct 
                        ie (fp+fp' in log space). Or in this case f[i][j] = fp[i] + fp[j] to minimize fc = fij + pi.


    Raises:
        null
    """

    # get the number of probes from the shape of the construct fitness matrix
    n = fc.shape[0]

    # subset the constuct fitness matrix to the upper triangle (its symmetric)
    # filter out bad constructs, based on w0 mask
    expressed = np.triu(w0).astype(np.bool)

    # initialize new weights matrix, assuming to start that all constructs are good
    w = np.ones((n,n))

    # filter weights
    w = filter_weights(w, w0)
    
    # get initial system of equations
    A, b = construct_system_ab(fc, w)
    
    # get initial solution
    fp, fij, eij = solve_iteration(A, b)

    # write status header
    if verbose:
        print("\t".join(["Iter","RMS","Relative Error"]))
    
    # init counter, relative error 
    counter = 0
    relative_error = 1
    while (relative_error > tol) & (counter < maxiter):
        # iterate the counter
        counter += 1
        
        # cache the current single probe fitnesses
        fp_old = fp.copy()

        # calculate normalized weights bases on pi scores
        w = biweight(eij, ag=ag, expressed=expressed)
        
        # filter weights
        w = filter_weights(w, w0)
        
        # get new weights and construct fitnesses
        A, b = construct_system_ab(fc, w)
        
        # get new solution
        fp, fij, eij = solve_iteration(A, b)

        # calculate relative error to the last iteration
        relative_error = np.sqrt( np.sum((fp - fp_old)**2) / np.max([1, sum(fp_old**2)]))
    
        # optionally print status to stdout 
        if verbose:
            j = np.sqrt(np.mean(fp**2))
            print("{}\t{:.4f}\t{:.6f}".format(counter,j,relative_error))

    # return final results
    return fp, fij, eij









if __name__ == "__main__":

    # set up argument parser
    parser = argparse.ArgumentParser(usage = globals()['__doc__'])
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose output")
    parser.add_argument("-f", "--construct_fitness", action="store", default=None, help="Path to construct fitness.")
    parser.add_argument("-w", "--construct_weights", action="store", default=None, help="Path to construct weights.")
    parser.add_argument("-o", "--output", action="store", default=None, \
                        help="Directory to write results to.")
    parser.add_argument("-a", "--n_stds", action="store", default=2, \
                        help="Number of standard deviation to use in Tukey biweight normalization")
    parser.add_argument("--tol", action="store", default=1e-3, \
                        help="Relative error tolerance")
    parser.add_argument("--maxiter", action="store", default=50, \
                        help="Maximum iterations to perform")
    # parse arguments
    options = parser.parse_args()

    if not options.construct_fitness and not options.construct_weights:
        raise BaseException("No input files provided.")

    # load input files
    fc = load_construct_fitnesses(options.construct_fitness)
    w0 = load_initial_weights(options.construct_weights)

    # compute
    fp, fij, eij = irls(fc, w0,
                        ag=options.n_stds,
                        tol=options.tol,
                        maxiter=options.maxiter,
                        verbose=options.verbose)

    if options.output:
        np.savetxt(os.path.join(options.output, "probe_fitnesses.csv"), fp, delimiter=",")
        np.savetxt(os.path.join(options.output, "pi_scores.csv"), eij, delimiter=",")
