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
import logging
import scipy.sparse as sps


def filter_weights(w, w0):
    """Filter weights for bad constructs and constructs with the only one target
    
    Set diagonal to zero and zero out specific entries
    Args:
        w (scipy sparse matrix): the matrix of weights to filter
        w0 (scipy sparse matrix): boolean matrix to zero out values from w
    Returns:
        w (scipy sparse matrix): the resulting filtered matrix

    """
    
    # replace the diagonal (ie same gene on both probes ) with zero ...
    # or do not include
    w.setdiag(0)

    # zero out previous determined bad constructs
    w = w.multiply(w0)
    return w

def biweight(eij, ag=2., expressed=None):
    """ Normalize weights according to a modified Tukey's Biweight
    see: http://mathworld.wolfram.com/TukeysBiweight.html
    Args:
        eij (scipy sparse matrix): a matrix of residuals from a fit
            of single probe fitnesses 
        ag (float): the number of standard deviation to use in normalization. 
            similar to (ag*sdtev)=c in tukey's biweight
        expressed (scipy sparse matrix): boolean matrix to filter residuals to
            in standard deviation calculation
    Returns:
        w (scipy sparse matrix): weights matrix
    """

    nonzero = eij.nonzero()
    l = nonzero[0].shape[0]
    ones = sps.csr_matrix((np.ones(l), nonzero), shape=eij.shape) 

    if expressed is None:
        # create a mask of all true
        expressed = ones.astype(np.bool)
        
    # calculate standard deviation of pi scores (residuals) for expressed
    # constructs
    sd = eij[expressed].std()
    
    # calculate new weights, aparently using a modified biweight model 
    # TODO, why not multiply by eij in front? 
    w = (ones - eij.power(2)/(ag*sd)**2).power(2)
    #w = (1-eij**2/(ag*sd)**2)**2
    
    # zero out anythin more than ag standard deviations from zero
    w[abs(eij)>(ag*sd)]=0

    return w

def construct_system_ab(fc, w, err=1e-6):
    """ Construct the system of equations to solve A*x=b.  Formulates A and b.

    Args:
        fc (scipy sparse matrix): matrix of consturct fitnesses
        w (scipy sprase matrix): matrix of weights
        err (float): small error term to apply to force diagonal not to be zero
    Returns:
        A (scipy sparse matrix)
        b (numpy array)
    """
    n = fc.shape[0]
    # initialize the imputed weights matrix
    A = w.copy()

    # replace diagnol with total number of constructs used for the given probe,
    # plus some small error term
    x = np.asarray(w.sum(axis=1)+err).reshape(-1)
    A.setdiag(np.asarray(w.sum(axis=1)+err).reshape(-1))

    # get new construct fitnesses to fit against initialize a vector with sum of
    # all construct fitnesses for every probe if the assumption of a normally
    # distributed genetic interactions about zero is met  this will be 0 or
    # close to it ... TODO: should it be the mean?

    b = np.array((fc.multiply(w)).sum(axis=1)).reshape((-1,))
    return A, b

def solve(A, b, exact=False):
    """ 
    Solve a system of equations corresponding to A*x=b.

    An underlying assumption of the dual knockout screens is that
    interations are rare. As such, the fitness of most double constructs will
    equal the multiplicative product of the single probe fitnesses. In log space
    this is an addative, fc = fa+fb. If a probe if measured against n other
    probes targeting different genetic elements then the sum of these constructs
    is expected to be 
    \sum_{i}^n fc_i = fc_1 + fc_2 + ... + fc_n = n*fa+fb_1+fb_2+fb3 

    As such the A matrix is constructed as integer weights and b is the sum of
    construct fitnesses for a given probe.
    ex:
    A = [ 3, 1, 1
          1, 3, 1,
          1, 1, 3]
    b = [ 
        ]

    Args:
        A (scipy sparse matrix): probe by probe matrix of weights 
        b (numpy vector): sum of constructs fitness by probe
        exact (bool): True if solving the system directly (square desgins),
            False otherwise, for example in multiple conditions
    Returns:
        x (numpy vector): vector of inputed probe fitness


    TODO: should we be using the fc matrix to get the nonzero instead of
    expressed? see issue #5
    """
    # convert to csr for computation
    #A = A.tocsr()

    # check to see if dimensions agree for solve
    r, c = A.shape
    if exact:
        if r!=c:
            raise RunTimeError("{}".format(
                "Trying to solve a non-square matrix. Useexact=False."))
        # find single probe fitnesses which satisfy expectation, you additive
        # property
        #x = np.linalg.solve(A, b)
        x = sps.linalg.spsolve(A, b)
    else:
        # if they do not agree then we must use an iterative least-squares
        # method
        x = sps.linalg.lsqr(A, b)[0]

    return x

def calc_eij(fp, fc, expressed=None):
    """ Calculate the residuals of a inputed probe based fitnesses.

    fij (matrix): A symmetric matrix of floats corresponding to the expected
    fitnesses for each contruct ie (fp+fp' in log space). Or in this case 
    f[i][j] = fp[i] + fp[j] to minimize fc = fij + pi.

    Args: 
        fp (numpy array): vector of inputed probe based fitnesses
        fc (scipy sparse matrix): probe by probe matrix of measured construct
            fitnesses
        expressed (scipy sparse matrix): boolean matrix mask indicating weather
            to include a construct in the calculation
    Returns:
        eij (scipy sparse matrix): the residuals of the expected constructs
        fitnesses bases on addative probe fitness from observed in the construct
        fitness matrix (fc)
    """
    # get indicies on nonzero elements
    if expressed is not None:
        nonzero  = expressed.nonzero()
    else:
        nonzero = fc.nonzero()

    # get expected double mutant fitness
    fij = sps.csr_matrix((fp[nonzero[0]] + fp[nonzero[1]], nonzero),
         shape=fc.shape)
    # get pi scores as observed deviations from expected
    eij = sps.triu(fc) - fij

    # make symmetric
    eij = eij + eij.transpose()

    return eij








def irls(fc, w0, ag=2, tol=1e-3, maxiter=50, all = True):
    """ The iterative least squares fitting function of single gene fitnesses
        
        Args:
            fc (matrix): A NxN matrix of observed fitnesses for each construct.
            w0 (matrix): A NxN matrix of boolean weights
            ag (float): Number of standard devs to use in biweight
            tol (float): The error tolerance threshold which stops the iteration.
            maxiter (int): max number of iteration to perform
        Return:
            guide_edit_distance (int) - number of one-nucleotide edits needed to transform the guide string into reference
    """ 
    # subset the constuct fitness matrix to the upper triangle (its symmetric)
    # filter out bad constructs, based on w0 mask
    expressed = sps.triu(w0).astype(np.bool)
   
    # write status header
    # if verbose:
    #     print("\t".join(["Iter","Relative Error"]))
    
    # init counter, relative error 
    counter = 0
    relative_error = 1
    while (relative_error > tol) & (counter < maxiter):
        if counter > 0: 
            # cache the current single probe fitnesses
            fp_old = fp.copy()

            # calculate normalized weights bases on pi scores
            w = biweight(eij, ag=ag, expressed=expressed)
        else:
            w = w0

        # filter weights
        w = filter_weights(w, w0)
        
        # get new weights and construct fitnesses
        A, b = construct_system_ab(fc, w)
        
        # get new solution
        fp = solve(A, b)
        eij = calc_eij(fp, fc, expressed)

        if counter > 0:
            # calculate relative error to the last iteration
            relative_error = np.sqrt( np.sum((fp - fp_old)**2) / np.max([1, sum(fp_old**2)]))
            logging.debug("Iteration: {}\t Relative Error: {:.6f}".format(counter,relative_error))

        # iterate the counter
        counter += 1

    # return final results
    return fp, eij




def irls_multi_condition(fc, w0, ag=2., tol=1e-3, maxiter=50, all = True):
    """ The iterative least squares fitting function of single gene fitnesses
        
        Args:
            fc (matrix): A NxN matrix of observed fitnesses for each construct.
            w0 (matrix): A NxN matrix of boolean weights
            ag (float): Number of standard devs to use in biweight
            tol (float): The error tolerance threshold which stops the iteration.
            maxiter (int): max number of iteration to perform
        Return:
            guide_edit_distance (int) - number of one-nucleotide edits needed to transform the guide string into reference
    """

    # init holders for cross condition system of equations
    Es = []
    # for each condition get the expressed constructs
    for i in range(len(fc)):
        # filter out bad constructs, based on w0 mask
        expressed = sps.triu(w0[i]).astype(np.bool)
        Es.append(expressed)

    # init counter, relative error 
    counter = 0
    relative_error = 1
    while (relative_error > tol) & (counter < maxiter):

        
        As, bs = [], []
        for i in range(len(fc)):
            
            if counter > 0:
                # calculate normalized weights bases on pi scores
                w = biweight(eijs[i], ag=ag, expressed=Es[i])
                # cache the current single probe fitnesses
                fp_old = fp.copy()
            else:
                w = w0[i]

            # filter weights
            w = filter_weights(w, w0[i])
             
            # get new weights and construct fitnesses
            A, b = construct_system_ab(fc[i], w)

            As.append(A)
            bs.append(b)

        # vstack the condition level systems
        A = sps.vstack(As)
        b = np.concatenate(bs)
    
        # get initial solution for the entire space
        fp = solve(A, b)

        eijs = []
        for i in range(len(fc)):
            eij = calc_eij(fp, fc[i], Es[i])
            eijs.append(eij)

        if counter > 0:
            # calculate relative error to the last iteration
            relative_error = np.sqrt( np.sum((fp - fp_old)**2) / np.max([1, sum(fp_old**2)]))
            logging.debug("Iteration: {}\t Relative Error: {:.6f}".format(counter,relative_error))

        # iterate the counter
        counter += 1

    # return final results
    return fp, eijs






