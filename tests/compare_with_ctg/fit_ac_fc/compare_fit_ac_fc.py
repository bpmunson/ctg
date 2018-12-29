import numpy as np
import pandas as pd

import os
import sys
sys.path.insert(0, "../../../ctg/core")

#from fit_ac_fc import *
from fit_ac_fc2 import Counts
from config import config

abundance_file = os.path.join(config.A549_test, "A549_abundance_thresholds.txt")
counts_file = os.path.join(config.A549_test, "A549_timepoint_counts.txt")
#times = np.array([[3,14, 21, 28], [3,14,21,28]])
times = np.array([3,14,21,28])

#ac, fc, allbad, sdfc, df, p_t, lfdr, names, lmbda, xfit, mask = fit_ac_fc(abundance_file, counts_file, times)
#ac, fc, allbad, sdfc, df, p_t, lfdr = fit_ac_fc(abundance_file, counts_file, times)
#ac, fc, allbad, sdfc, p_t, names = fit_ac_fc(abundance_file, counts_file, times)

filesList = ['a1_ctg.csv',
             'a2_ctg.csv',
             'allbad_ctg.csv',
             'df_ctg.csv',
             'fc_ctg.csv',
             'lfdr_fc_ctg.csv',
             'p_t_ctg.csv',
             'sdfc_ctg.csv',
             'lambda1_ctg.csv',
             'lambda2_ctg.csv',
             'xfit1_ctg.csv',
             'xfit2_ctg.csv',
             'good1_ctg.csv',
             'good2_ctg.csv']

# for fn in filesList:
#     pd.read_csv(fn)

def compare_a():
    a1_df = pd.read_csv(filesList[0])
    a2_df = pd.read_csv(filesList[1])

    assert (a1_df.iloc[:,0] != a2_df.iloc[:,0]).sum() == 0

    ac_ctg = pd.DataFrame([a1_df.iloc[:,1], a2_df.iloc[:,1]],).T
    ac_ctg.index = a1_df.iloc[:,0]

    assert np.allclose(ac_ctg.as_matrix(), ac.T)

def comapre_fc():
    fc_ctg = pd.read_csv(filesList[4], index_col=0)

    assert np.allclose(fc_ctg.as_matrix().ravel(), fc)

def compare_allbad():
    allbad_ctg = pd.read_csv(filesList[2], index_col=0)

    assert (allbad_ctg.as_matrix().ravel() == np.array(allbad)).all()

def compare_sdfc():
    sdfc_ctg = pd.read_csv(filesList[7], index_col=0)

    # print(sdfc)
    # print(sdfc_ctg.as_matrix().ravel())

    assert np.allclose(sdfc, sdfc_ctg.as_matrix().ravel())

def compare_df():
    df_ctg = pd.read_csv(filesList[3], index_col=0)

    assert np.allclose(df_ctg.as_matrix().ravel(), df)

def compare_p_t():
    p_t_ctg = pd.read_csv(filesList[6], index_col=0)

    assert np.allclose(p_t_ctg.as_matrix().ravel(), p_t)

def compare_lfdr():
    lfdr_ctg = pd.read_csv(filesList[5], index_col=0)

    print(lfdr_ctg.as_matrix().ravel())
    print(lfdr)

    print(np.allclose(lfdr_ctg.as_matrix().ravel(), lfdr))

def _compare_lambda():
    #Note: lambda is usually not returned. Use for debugging
    lmbda1_df = pd.read_csv(filesList[8], index_col=0)
    lmbda2_df = pd.read_csv(filesList[9], index_col=0)

    lmbda_ctg = pd.DataFrame([lmbda1_df.iloc[:,0], lmbda2_df.iloc[:,0]],).T

    assert np.allclose(lmbda_ctg.as_matrix(), lmbda.T)

def _compare_xfit():
    xfit1_df = pd.read_csv(filesList[10], index_col=0)
    xfit2_df = pd.read_csv(filesList[11], index_col=0)

    assert np.allclose(xfit1_df.as_matrix(), xfit[0])
    assert np.allclose(xfit2_df.as_matrix(), xfit[1])

def _compare_mask():
    good1_df = pd.read_csv(filesList[12], index_col=0)
    good2_df = pd.read_csv(filesList[13], index_col=0)

    assert np.allclose(good1_df.as_matrix(), mask[0])
    assert np.allclose(good2_df.as_matrix(), mask[1])

if __name__ == "__main__":
    def test_suite(): 
        compare_a()
        comapre_fc()
        compare_allbad()
        compare_sdfc()
        #compare_df()
        #compare_p_t()
        #compare_lfdr()

        # _compare_lambda()
        # _compare_xfit()
        # _compare_mask()

    def fit_ac_fc(abundance_file, counts_file, times): 
        if isinstance(counts_file, str): 
            c = Counts.from_file(counts_file)
            c.fit_ac_fc(pd.read_csv(abundance_file, sep='\t', index_col=0))

        else: 
            c = Counts(counts_file)
            c.fit_ac_fc(pd.read_csv(abundance_file, sep='\t', index_col=0))

        return c.ac, c.fitness, c.mask.allbad, c.sdfc, c.p_t, c.names

    global ac, fc, allbad, sdfc, p_t, names

    #Testing files input
    ac, fc, allbad, sdfc, p_t, names = fit_ac_fc(abundance_file, counts_file, times)
    test_suite()

    print('Passed file inputs (implicit)')

    # #Testing files input
    # ac, fc, allbad, sdfc, p_t, names = fit_ac_fc(abundance_file, counts_file, times, method='explicit', 
    #                                                 n_reps=2, 
    #                                                 columns_map=[[0,1],[2,3],[4,5],[6,7]])
    # test_suite()    
    # print('Paseed file inputs (explicit)')

    #Testing dataframe input
    counts_df = pd.read_csv(counts_file, sep='\s+')
    ac, fc, allbad, sdfc, p_t, names = fit_ac_fc(abundance_file, counts_df, times)
    test_suite()

    print('Passed dataframe input (implicit)')

    #Testing explicit method
    # ac, fc, allbad, sdfc, p_t, names = fit_ac_fc(abundance_file, counts_df, times, method='explicit', 
    #                                                 n_reps=2, 
    #                                                 columns_map=[[0,1],[2,3],[4,5],[6,7]])
    # test_suite()
    # print('Passed dataframe input (explicit')

    #TODO: Test changing the column order
