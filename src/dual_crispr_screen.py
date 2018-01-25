
import os
import argparse
import pandas as pd
import numpy as np

import fit_ac_fc
import irls
import weighted_pi
import bootstrap_pi


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


class Screen(object):

	def __init__(self, abundance_file, timepoint_counts_file, times):
		self.abundance_file = abundance_file
		self.timepoint_counts_file = timepoint_counts_file
		self.times = times

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
		self.fij = None
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
																		 keep_names=True)
		

		# store results
		self.names = names
		self.fc = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], fc, return_dataframe=True)
		self.allbad = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], allbad, return_dataframe=True)
		self.sdfc = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], sdfc, return_dataframe=True)
		self.pp = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], p_t, return_dataframe=True)

		# get initial weights
		w0 = self._get_initial_weights(allbad)
		self.w0 = self._make_symmetric_matrix(names['probe_a_id'], names['probe_b_id'], w0, return_dataframe=True).astype(int)

		# store probes
		self.probes = np.array(self.fc.index)

		# get and store target ids
		self.targets = self._build_target_array()

		if self.options.testing:
			# reorder according to RS/AB for testing purposes
			# this is required because of the random number generator needs the same format
			benchmark = '~/crappy/data/test_data/output_data/Notebook8Test_fc_0_benchmark.csv'
			benchmark = pd.read_csv(benchmark, sep=",", header=0, index_col=0)

			# reset indices 
			self.fc = self.fc.reindex(index=benchmark.index, columns=benchmark.index)
			self.allbad = self.allbad.reindex(index=benchmark.index, columns=benchmark.index)
			self.sdfc = self.sdfc.reindex(index=benchmark.index, columns=benchmark.index)
			self.pp = self.pp.reindex(index=benchmark.index, columns=benchmark.index)
			self.w0 = self.w0.reindex(index=benchmark.index, columns=benchmark.index)
			# TODO
			# also need to reload the posteriors because they are off from fit_ac_fc at the moment
			path = '../data/test_data/output_data/Notebook8Test_pp_0_benchmark.csv'
			self.pp = pd.read_csv(path, sep=",", index_col=0, header=0)

			# restore probes
			self.probes = list(self.fc.index)
			self.targets = self._build_target_array()

	def run_pi_score_calculation(self):
		""" Description
		"""
		if self.fc is None or self.w0 is None:
			raise AssertionError('No construct fitness or weights available. Must first need to run construct fitting.')
		
		# run irls
		fp, fij, eij = irls.irls(self.fc.values,
								 self.w0.values,
								 ag = self.options.n_stds,
								 tol = self.options.tol,
								 maxiter = self.options.maxiter,
								 verbose = self.options.verbose)
		
		# store results
		self.fp = fp
		self.fij = fij
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
													self.w0.values,
													self.probes,
													self.targets,
													n_probes_per_target=self.options.n_probes_per_target,
													null_target_id = self.options.null_target_id,
													null = self.options.null_aware,
													fp_0 = None
													)

		self.pi_scores = pi_scores
		self.target_fitness = target_fitness

	def run_bootstrap(self):
		""" Descriptions
		"""
		pi_iter, fitness_iter = bootstrap_pi.run_iteration(
															self.fc.values,
															self.pp.values,
															self.sdfc.values,
															self.w0.values,
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
															use_full_dataset_for_ranking = self.options.use_full_dataset_for_ranking	
														  )

		# store results 
		self.pi_scores_iter = pi_iter
		self.fitness_iter = fitness_iter


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

		# get a list of all probes to targets from names datafame
		a = self.names[['probe_a_id','target_a_id']]
		b = self.names[['probe_b_id','target_b_id']]
		a.columns = ['probe_id','target_id']
		b.columns = ['probe_id','target_id']
		c = pd.concat([a,b],axis=0).drop_duplicates()

		# merge together and take target ids
		targets = np.array(pd.merge(pd.DataFrame(self.probes,columns=['probe_id']), c)['target_id'])

		return targets