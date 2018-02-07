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
import scipy.sparse as sps


def load_matrix_csv(fp, sep=",", index_col=0, header=0):
    """ Loads precomputed constuct fitnesses from a csv file
    """

    fc = pd.read_csv(fp, sep=sep, header=header, index_col=index_col).values
    return fc

def filter_weights(w, w0):
    """Filter weights for bad constructs and constructs with the only one target
    """
    
    # replace the diagonal (ie same gene on both probes ) with zero ... or do not include
    #np.fill_diagonal(w, 0)
    w.setdiag(0)

    # zero out previous determined bad constructs 
    w = w.multiply(w0)

    return w

def biweight(eij, ag=2, expressed=None):
    """ Normalize weights according to a modified Tukey's Biweight
    """

    nonzero = eij.nonzero()
    l = nonzero[0].shape[0]
    ones = sps.csr_matrix((np.ones(l), nonzero), shape=eij.shape) 

    if expressed is None:
        # create a mask of all true
        expressed = ones.astype(np.bool)
        
    # calculate standard deviation of pi scores (residuals) for expressed constructs
    sd = eij[expressed].std()
    
    # calculate new weights, aparently using a modified biweight model 
    # TODO, why not multiply by eij in front? 
    w = (ones - eij.power(2)/(ag*sd)**2).power(2)
    #w = (1-eij**2/(ag*sd)**2)**2
    
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
    A.setdiag(np.asarray(w.sum(axis=1)+err).reshape(-1))

    # get new construct fitnesses to fit against
    # initialize a vector with sum of all construct fitnesses for every probe
    # if the assumption of a normally distributed genetic interactions about zero is met 
    # this will be 0 or close to it ... TODO: should it be the mean? 
    b = np.array((fc.multiply(w)).sum(axis=1)).reshape((-1,))
    return A, b

def solve_iteration(A, b, fc, expressed, all=True):
    """ 
    Solves for the individual probe fitnesses

    TODO: should we be using the fc matrix to get the nonzero instead of expressed? see issue #5
    """
    # convert to csr for computation
    A = A.tocsr()

    # check to see if dimensions aggree for solve
    r, c = A.shape
    if r==c:
        # find single probe fitnesses which satisfy expectation, you additive property
        #x = np.linalg.solve(A, b)
        x = sps.linalg.spsolve(A, b)
    else:
        # if they do not agree then we must use an iterative least-squares method
        x = sps.linalg.lsqr(A, b)[0]

    # TODO: by taking 
    # get indicies on nonzero elements
    if all:
        nonzero  = expressed.nonzero()
    else:
        nonzero = fc.nonzero()

    # get expected double mutant fitness
    fij = sps.csr_matrix((x[nonzero[0]] + x[nonzero[1]], nonzero), shape=fc.shape)
    # get pi scores as observed deviations from expected
    eij = sps.triu(fc) - fij

    # make symmetric
    fij = fij + fij.transpose()
    eij = eij + eij.transpose()

    return x, fij, eij



def irls(fc, w0, ag=2, tol=1e-3, maxiter=50, verbose=False, all = True):
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

    # subset the constuct fitness matrix to the upper triangle (its symmetric)
    # filter out bad constructs, based on w0 mask
    #expressed = np.triu(w0).astype(np.bool)
    expressed = sps.triu(w0).astype(np.bool)
    # initialize new weights matrix, assuming to start that all constructs are good
    #w = np.ones((n,n))

    # filter weights
    w = filter_weights(w0, w0)
    
    # get initial system of equations
    A, b = construct_system_ab(fc, w)
    
    # get initial solution
    fp, fij, eij = solve_iteration(A, b, fc, expressed, all=all)

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
        fp, fij, eij = solve_iteration(A, b, fc, expressed, all=all)

        # calculate relative error to the last iteration
        relative_error = np.sqrt( np.sum((fp - fp_old)**2) / np.max([1, sum(fp_old**2)]))
    
        # optionally print status to stdout 
        if verbose:
            j = np.sqrt(np.mean(fp**2))
            print("{}\t{:.4f}\t{:.6f}".format(counter,j,relative_error))

    # return final results
    return fp, fij, eij




