import os
import argparse
import pandas as pd
import numpy as np
import fit_ac_fc
import irls
import weighted_pi
import bootstrap_pi
import pickle
import scipy.sparse as sps


class Options(object):
    def __init__(self):

        ###########################################
        # defaults
        ###########################################
        # fit_ac_fc
        self.n_good = 2
        self.replicate_axis = 0
        self.samples_axis = 1
        self.timepoints_axis = 2
        # irls
        self.verbose = False
        self.n_stds = 2
        self.tol = 1e-3
        self.maxiter = 50
        # weighted pi
        self.n_probes_per_target = 2
        self.epsilon = 1e-6
        self.null_target_id = "0"
        self.null_aware = True
        self.use_full_dataset_for_ranking = True
        self.niter = 2
        self.testing = False
        self.output = None
        self.all = True




class Screen(object):

    def __init__(self, timepoint_counts_file, times,
        abundance_file=None, min_counts_threshold=10, verbose=False):

        self.abundance_file = abundance_file
        self.timepoint_counts_file = timepoint_counts_file
        self.times = times
        self.min_counts_threshold=min_counts_threshold
        self.verbose=verbose

        # initialize default options
        # TODO: allow for argparse, config, and setting of options upon initialization
        self.options = Options()

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
        # parse the times argument if provide
        ac, fc, allbad, sdfc, df, p_t, lfdr, names = fit_ac_fc.fit_ac_fc(self.abundance_file,
                                                                         self.timepoint_counts_file,
                                                                         self.times,
                                                                         n_good = self.options.n_good,
                                                                         replicate_axis = self.options.replicate_axis,
                                                                         samples_axis = self.options.samples_axis,
                                                                         timepoints_axis = self.options.timepoints_axis,
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

        if self.options.testing:
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

    def make_single_gene_screen(self):
        nontargeting = [i for i in list(self.fc.columns) if i.startswith("NonTargeting")]
        everything_else = [i for i in list(self.fc.columns) if not i.startswith("NonTargeting")]

        self.fc[everything_else] = 0
        self.pp[everything_else] = 0
        self.sdfc[everything_else] = 0
        self.w0[everything_else] = 0
        self.fc[everything_else] = 0

    def run_pi_score_calculation(self):
        """ Description
        """
        if self.fc is None or self.w0 is None:
            raise AssertionError('No construct fitness or weights available. Must first need to run construct fitting.')
        
        # run irls
        fp, eij = irls.irls(self.fc,
                                 self.w0,
                                 ag = self.options.n_stds,
                                 tol = self.options.tol,
                                 maxiter = self.options.maxiter,
                                 verbose = self.options.verbose)
        

        # store results
        self.fp = fp
        self.eij = eij

    def run_weighting(self):
        """ Descriptions
        """
        if self.fc is None or self.w0 is None:
            raise AssertionError('No construct fitness or weights available. Must first need to run construct fitting.')
        # compute weighted estimates
        pi_scores, target_fitness = weighted_pi.weight_by_target(
                                                    self.eij, 
                                                    self.fp,
                                                    self.w0,
                                                    self.probes,
                                                    self.targets,
                                                    n_probes_per_target=self.options.n_probes_per_target,
                                                    null_target_id = self.options.null_target_id,
                                                    null = self.options.null_aware,
                                                    pre_computed_ranks = None
                                                    )

        self.pi_scores = pi_scores
        self.target_fitness = target_fitness

    def run_bootstrap(self):
        """ Descriptions
        """
        pi_iter, fitness_iter = bootstrap_pi.run_iteration(
                                                            self.fc,
                                                            self.pp,
                                                            self.sdfc,
                                                            self.w0,
                                                            self.probes,
                                                            self.targets,
                                                            ag = self.options.n_stds,
                                                            tol = self.options.tol,
                                                            maxiter = self.options.maxiter,
                                                            n_probes_per_target = self.options.n_probes_per_target,
                                                            epsilon = self.options.epsilon,
                                                            null_target_id = self.options.null_target_id,
                                                            null = self.options.null_aware,
                                                            niter = self.options.niter,
                                                            verbose = self.options.verbose,
                                                            testing = self.options.testing,
                                                            use_full_dataset_for_ranking = self.options.use_full_dataset_for_ranking,
                                                            all=self.options.all    
                                                          )

        # store results 
        self.pi_scores_iter = pi_iter
        self.fitness_iter = fitness_iter
        self.results = bootstrap_pi.summarize(self.pi_scores_iter, self.fitness_iter)

    def pickle(self):
        output = self.options.output
        pickle.dump(self, open(output, 'wb'))

    def summarize(self):
        """ Description
        """
        results = bootstrap_pi.summarize(self.pi_scores_iter, self.fitness_iter)

        # store results
        self.results = results



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


