import os
import numpy as np
import pandas as pd
import ctg
from ctg.core.config import config
import logging

def _compare_to_benchmark(testing, reference):
    """ Compare a numpy array to the benchmark """
    assert np.allclose(testing, reference)

def _compare_probe_fitnesses(scorer, file):
    # compare probe fitnesses
    reference = pd.read_csv(file, sep="\t", header=0)['fitness'].values
    _compare_to_benchmark(scorer.fp, reference)

def _compare_pi_scores(scorer, file):
    # compare to pi scores 
    reference = pd.read_csv(file, sep="\t",header=0, index_col=0).values
    _compare_to_benchmark(scorer.eij.todense(), reference)

def _compare_weighted_pi_scores(scorer, file):
    # compare weighted pi scores
    reference = pd.read_csv(file, sep="\t", header=0)['pi'].values
    _compare_to_benchmark(scorer.pi_scores['pi'].values, reference)

def _compare_weighted_fitnesses(scorer, file):
    # compare weighted target fitnesses
    reference = pd.read_csv(file, sep="\t", header=0)['fitness'].values
    _compare_to_benchmark(scorer.target_fitness[0].values, reference)

def _compare_pi_scorer_iter(scorer, file):
    # compare iter pi_scores 
    reference = pd.read_csv(file, sep="\t", header=0, index_col=[0,1]).values
    _compare_to_benchmark(scorer.pi_scores_iter.values, reference)

def _compare_fitness_iter(scorer, file):
    # compare fitness iter
    reference = pd.read_csv(file, sep="\t", header=0, index_col=0).values
    _compare_to_benchmark(scorer.fitness_iter.values, reference)

def _compare_results(scorer, file):
    # compare final results file
    reference = pd.read_csv(file, sep="\t", header=0, index_col=[0,1]).values
    _compare_to_benchmark(scorer.results.iloc[:,2:].values, reference)

def main():
    """ Run calculation and perform tests"""
    files = {
        "fp": "benchmarks/testing_fp_benchmark.csv",
        "eij": "benchmarks/testing_eij_benchmark.csv",
        "weighted_pi_scores": "benchmarks/testing_weighted_pi_scores_benchmark.csv",
        "weighted_fitness": "benchmarks/testing_weighted_finesses_benchmark.csv",
        "pi_scores_iter": "benchmarks/testing_pi_scores_iter_benchmark.csv",
        "fitness_iter": "benchmarks/testing_fitness_iter_benchmark.csv",
        "results": "benchmarks/testing_results_benchmark.csv",
    }   

    abundance_file = os.path.join(config.A549_test, "A549_abundance_thresholds.txt")
    counts_file = os.path.join(config.A549_test, "A549_timepoint_counts.txt")
    times = '3,14,21,28'

    logging.info("Initializing Scorer")
    # initialize scorer
    scorer = ctg.Scorer(
        counts_file,
        times,
        abundance_file= abundance_file,
        niter = 2,  
        null_target_id = "NonTargetingControl",
        testing = True)

    # run calculation with current build
    scorer.run_construct_fitting()
    scorer.run_pi_score_calculation()
    scorer.run_weighting()
    scorer.run_sampling()
    scorer.summarize()

    logging.info("Running Tests.")

    # test results
    _compare_probe_fitnesses(scorer, files['fp'])
    _compare_pi_scores(scorer, files['eij'])
    _compare_weighted_pi_scores(scorer, files['weighted_pi_scores'])
    _compare_weighted_fitnesses(scorer, files['weighted_fitness'])
    _compare_pi_scorer_iter(scorer, files['pi_scores_iter'])
    _compare_fitness_iter(scorer, files['fitness_iter'])
    _compare_results(scorer, files['results'])

    return 0

if __name__ == "__main__":
    # set up logger
    logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s',
                        level=logging.INFO) 
    # Run Tests
    ret_code = main()

    if ret_code == 0:
        logging.info("All Tests Passed")
