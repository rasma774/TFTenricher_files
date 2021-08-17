
import pandas as pd
import numpy as np
import os
import scipy.stats as sts

import TFTenricher

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2021, Sävsjö'
__contact__ = 'rasma774@gmail.com'


def benjaminihochberg_correction(p, FDR=0.05):
    """


    Parameters
    ----------
    p : np.array
        array of p-values of independent tests.
    FDR : float, optional
        False discovey rate. 0 < FDR < 1. The default is 0.05.

    Returns
    -------
    passes_FDR
        A vector of same length as input parameter p, with bool value True
        where the test passed a BH FDR correction.

    """
    sorted_p = np.sort(p)
    rank = np.arange(1, len(p)+1)

    BH_crit = (rank/rank[-1])*FDR

    if np.any(sorted_p < BH_crit):
        thresh = sorted_p[sorted_p < BH_crit][-1]
        return p <= thresh
    else:
        return p < 0

TFnames = pd.read_csv('tfnames.txt', header=None)[0].values
diseases = os.listdir('expression_atlas')

methods = ['GO', 'KEGG', 'REACTOME', 'GWAS']
res_target = {}
res_TF = {}
for d in diseases:
    pvals_tmp = pd.read_csv('expression_atlas/' + d + '/TFp.csv').iloc[:,-1].values

    is_sign = benjaminihochberg_correction(pvals_tmp)
    if (is_sign.sum() > 300) or (is_sign.sum() < 10):
        continue
    TFtmp = TFnames[is_sign]

    enr = TFTenricher.TFTenricher(TFtmp, mapmethod='corr')

    restemp_target = {}
    restemp_TF = {}
    for method_name in methods:
        enr.downstream_enrich(db=method_name)
        enr.enrichments.to_csv('results/' + method_name + '/' + d + '.csv')
        restemp_target[method_name] = enr.enrichments.FDR.sum()


    enr.target_genes = TFtmp
    for method_name in methods:
        enr.downstream_enrich(db=method_name)
        enr.enrichments.to_csv('results/' + method_name + '/' + d + '_TFs.csv')
        restemp_TF[method_name] = enr.enrichments.FDR.sum()

    res_target[d] = restemp_target
    res_TF[d] = restemp_TF

res_target = pd.DataFrame(res_target)
res_target.to_csv('disease_target_res_after_review.csv')

res_TF = pd.DataFrame(res_TF)
res_TF.to_csv('disease_TF_res_after_review.csv')

print(res_target.index)
for i in range(4):
    print(sts.wilcoxon(res_target.values[i,:], res_TF.values[i,:], alternative='greater'))
print(np.median(res_target.values, 1))
print(np.median(res_TF.values, 1))

print(np.mean(res_target.values, 1))
print(np.mean(res_TF.values, 1))
