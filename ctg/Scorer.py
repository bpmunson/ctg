import os
import argparse
import logging
import pickle

import pandas as pd
import numpy as np
import scipy.sparse as sps

import ctg.fit_ac_fc as fit_ac_fc
import ctg.irls as irls
import ctg.weight as weight
import ctg.sample as sample

class Scorer(object):

    def __init__(self, timepoint_counts_file, times,
        abundance_file=None,
        min_counts_threshold=10,
        min_time_points = 2,
        replicate_axis = 0,
        samples_axis = 1,
        timepoints_axis = 2,
        verbose = False,
        bi_weight_stds = 2,
        tol = 1e-3,
        maxiter = 50,
        n_probes_per_target = 2,
        epsilon = 1e-6,
        null_target_id = "NonTargeting",
        use_full_dataset_for_ranking = True,
        niter = 2,
        testing = False,
        output = None,
        pickle_output = None):


        if isinstance(times, str):
            times = np.array([int(i) for i in times.split(",")])

        self.timepoint_counts_file = timepoint_counts_file
        self.times = times
        self.abundance_file = abundance_file
        self.min_counts_threshold = min_counts_threshold 
        self.min_time_points = min_time_points
        self.replicate_axis = replicate_axis
        self.samples_axis = samples_axis
        self.timepoints_axis = timepoints_axis
        self.verbose = verbose
        self.bi_weight_stds = bi_weight_stds
        self.tol = tol
        self.maxiter = maxiter
        self.n_probes_per_target = n_probes_per_target
        self.epsilon = epsilon
        self.null_target_id = null_target_id
        #self.use_full_dataset_for_ranking = use_full_dataset_for_ranking
        self.niter = niter
        self.testing = testing
        self.output = output
        self.pickle_output = pickle_output

        ###########################################
        # results
        ###########################################
        self.probes = None
        self.targets = None
        self.construct_fitnesses = None
        self.names = None
        self.fc = None
        self.allbad = None
        self.sdfc = None
        self.pp = None
        self.fp = None
        self.eij = None
        self.pi_scores = None
        self.target_fitness = None
        self.pi_scores_iter = None
        self.results = None

    def __repr__(self):
        return '{}: {}'.format(self.__class__.__name__, self.abundance_file)

    def run_construct_fitting(self):
        """ Description
        """
        logging.info("Calculating construct fitnesses.")
        # parse the times argument if provide
        ac, fc, allbad, sdfc, df, p_t, lfdr, names = fit_ac_fc.fit_ac_fc(self.abundance_file,
                                                                         self.timepoint_counts_file,
                                                                         self.times,
                                                                         n_good = self.min_time_points,
                                                                         replicate_axis = self.replicate_axis,
                                                                         samples_axis = self.samples_axis,
                                                                         timepoints_axis = self.timepoints_axis,
                                                                         keep_names=True,
                                                                         min_counts_threshold=self.min_counts_threshold,
                                                                         verbose=self.verbose)


        self.fc0 = fc
        # store results
        self.names = names
         # get initial weights
        w0 = self._get_initial_weights(allbad)

        # build into sparse matrices
        self.fc, self.probes = self._make_sparse_matrix(names['probe_a_id'], names['probe_b_id'], fc, return_probes=True)
        self.sdfc = self._make_sparse_matrix(names['probe_a_id'], names['probe_b_id'], sdfc)
        self.pp = self._make_sparse_matrix(names['probe_a_id'], names['probe_b_id'], p_t)
        self.w0 = self._make_sparse_matrix(names['probe_a_id'], names['probe_b_id'], w0)
        # get and store target ids
        self.targets = self._build_target_array()

        if self.testing:
            # get dataframes
            self.fc = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], fc, return_dataframe=True)
            self.sdfc = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], sdfc, return_dataframe=True)
            self.pp = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], p_t, return_dataframe=True)
            self.w0 = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], w0, return_dataframe=True).astype(int)

            # also need to reload the posteriors because they are off from fit_ac_fc at the moment
            src = os.path.dirname(os.path.realpath(__file__))
            path = os.path.join(src,"..",'data','test_data','output_data','Notebook8Test_pp_0_benchmark.csv')
            self.pp = pd.read_csv(path, sep=",", index_col=0, header=0)
            
            # reorder according to RS/AB for testing purposes
            # this is required because of the random number generator needs the same format
            benchmark = self.pp

            # reset indices
            self.fc = self.fc.reindex(index=benchmark.index, columns=benchmark.index)
            self.sdfc = self.sdfc.reindex(index=benchmark.index, columns=benchmark.index)
            self.pp = self.pp.reindex(index=benchmark.index, columns=benchmark.index)
            self.w0 = self.w0.reindex(index=benchmark.index, columns=benchmark.index)

            # restore probes
            self.probes = list(self.fc.index)
            self.targets = self._build_target_array()

            # convert back to sparse matrix
            self.fc = sps.csr_matrix(self.fc)
            self.pp = sps.csr_matrix(self.pp)
            self.w0 = sps.csr_matrix(self.w0)
            self.sdfc = sps.csr_matrix(self.sdfc)

    def run_pi_score_calculation(self):
        """ Description
        """
        logging.info("Calculating single probe fitnesses and probe level interaction scores.")
        if self.fc is None or self.w0 is None:
            raise AssertionError('No construct fitness or weights available. Must first need to run construct fitting.')
        
        # run irls
        fp, eij = irls.irls(self.fc,
                                 self.w0,
                                 ag = self.bi_weight_stds,
                                 tol = self.tol,
                                 maxiter = self.maxiter)
        

        # store results
        self.fp = fp
        self.eij = eij

    def run_weighting(self):
        """ Descriptions
        """
        logging.info("Weighting probe level interactions.")
        if self.fc is None or self.w0 is None:
            raise AssertionError('No construct fitness or weights available. Must first need to run construct fitting.')
        # compute weighted estimates
        pi_scores, target_fitness = weight.weight_by_target(
                                                    self.eij, 
                                                    self.fp,
                                                    self.w0,
                                                    self.probes,
                                                    self.targets,
                                                    n_probes_per_target=self.n_probes_per_target,
                                                    null_target_id = self.null_target_id,
                                                    pre_computed_ranks = None
                                                    )

        self.pi_scores = pi_scores
        self.target_fitness = target_fitness

    def run_sampling(self):
        """ Descriptions
        """
        logging.info("Calculating interaction scores by iterative subsampling.")
        pi_iter, fitness_iter = sample.run_iteration(
                                                            self.fc,
                                                            self.pp,
                                                            self.sdfc,
                                                            self.w0,
                                                            self.probes,
                                                            self.targets,
                                                            ag = self.bi_weight_stds,
                                                            tol = self.tol,
                                                            maxiter = self.maxiter,
                                                            n_probes_per_target = self.n_probes_per_target,
                                                            epsilon = self.epsilon,
                                                            null_target_id = self.null_target_id,
                                                            niter = self.niter,
                                                            testing = self.testing,
                                                            )

        # store results 
        self.pi_scores_iter = pi_iter
        self.fitness_iter = fitness_iter
        self.results = sample.summarize(self.pi_scores_iter, self.fitness_iter)

    def pickle(self):
        output = self.pickle_output
        print(output)
        pickle.dump(self, open(output, 'wb'))

    def summarize(self):
        """ Description
        """
        logging.info("Summarizing results.")
        results = sample.summarize(self.pi_scores_iter, self.fitness_iter)

        # store results
        self.results = results

        if self.output:
            self.results.to_csv(self.output, sep="\t", header=True, index=False)





    ###############################################################################################
    ### Helper functions
    ###############################################################################################
    def _make_symmetric_matrix(self, probe_a, probe_b, feature, return_dataframe=False):
        """ make a symmetric matrix or data frame from a list of probes pairs and feature vector
        """
        #TODO: assert they are they same length

        # make a dataframe
        df = pd.concat([probe_a, probe_b], axis=1)
        df.loc[:,'feature'] = feature
        df.columns=['probe_a_id', 'probe_b_id', 'feature']

        # pivot the table
        dfp = df.pivot(index="probe_a_id", columns="probe_b_id", values='feature')

        # get a common index
        index = dfp.index.union(dfp.columns)

        # reindex the pivoted table
        dfp = dfp.reindex(index=index, columns=index, fill_value = 0).fillna(0)

        # make symmetric
        vals = dfp.values
        vals = vals + np.transpose(vals)

        if return_dataframe:
            return pd.DataFrame(vals, index=dfp.index, columns=dfp.columns)
        else:
            return vals

    def _make_sparse_matrix(self, probe_a, probe_b, feature, return_probes=False):
        """ The iterative least squares fitting function of single gene fitnesses

        Args:
            probe_a (list): A list of a probes corresponding to one of the crispr guides in a dual knockout
            probe_a (list): A list of a probes corresponding to the other of the crispr guides in a dual knockout
            feature (array): An array of features corresponding to the probe pair defined by probe_a and probe_b
        Returns:
            m (scipy sparse matrix): a sparse matrix with all probes by all probes, populated with the feature vector
            union (list): a list of all the probes in the screen corresponding to the sparse matrix
        Raises:
        """
        # make union of probe pairs
        #union = sorted(set(probe_a + probe_b))
        union = np.unique(np.concatenate((probe_a, probe_b)))
        # enumerate and hash
        union_labels = {j:i for i,j in enumerate(union)}
        # get total probes in screen
        l = len(union)
        # get integer labels for the probe lists
        probe_a_i = [union_labels[i] for i in probe_a]
        probe_b_i = [union_labels[j] for j in probe_b]
        # construct sparse matrix
        m = sps.csr_matrix((feature, (probe_a_i, probe_b_i)), shape=(l,l))
        # make symmetric
        m = m + m.transpose()
        # return sparse matrix and probes list
        if return_probes:
            return m, union
        else:
            return m

    def _get_initial_weights(self, allbad):

        """ Make intitial boolean weights from weather all replicates were bad for a given construct
        """
        if allbad is None:
            raise AssertionError('Cannot get initial weights without first running construct fitting.')

        # convert {0,1} allbad labels to boolean and take the opposite
        res = (~allbad).astype(int)

        return res

    def _build_target_array(self):
        """ Description
        """

        if self.probes is None:
            raise AssertionError('Must first construct probe array.')
        # get a list of all probes to targets from names datafame
        a = self.names[['probe_a_id','target_a_id']]
        b = self.names[['probe_b_id','target_b_id']]
        a.columns = ['probe_id','target_id']
        b.columns = ['probe_id','target_id']
        c = pd.concat([a,b],axis=0).drop_duplicates()

        # merge together and take target ids
        targets = np.array(pd.merge(pd.DataFrame(self.probes,columns=['probe_id']), c)['target_id'])

        return targets

    def _check_symmetric(a, tol=1e-8):
        return np.allclose(a, a.T, atol=tol,equal_nan=True)


