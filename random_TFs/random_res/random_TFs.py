import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import TFTenricher

import sys
sys.path.append('../src/')

import map2trgt_utils

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2021, Sävsjö'
__contact__ = 'rasma774@gmail.com'

all_tfs = pd.read_csv('tfnames.txt', header=None).values.T[0]

np.random.seed(0)



for mapmethod in ('corr', 'TRRUST', 'PPI'):
    n_FDR = []

    if mapmethod == 'corr':
        mapmethod_name = 'corr'
    elif mapmethod == 'TRRUST':
        mapmethod_name=map2trgt_utils.trrust_genes
    elif mapmethod == 'PPI':
        mapmethod_name=map2trgt_utils.STRING_ppi

    for i in range(100):
        print(i)
        pws_tmp = []
        TF_rand = np.random.choice(all_tfs, size=100, replace=False)
        enr = TFTenricher.TFTenricher(TF_rand, mapmethod=mapmethod_name, silent=True, top_n_genes=200)
        enr.downstream_enrich()
        corr_false_GO = enr.enrichments.FDR.sum()

        enr.downstream_enrich(db='REACTOME')
        corr_false_REACTOME = enr.enrichments.FDR.sum()

        enr.downstream_enrich(db='KEGG')
        corr_false_KEGG = enr.enrichments.FDR.sum()

        enr.downstream_enrich(db='GWAS')
        corr_false_GWAS = enr.enrichments.FDR.sum()

        n_FDR.append(
            {
                'GO':corr_false_GO,
                'KEGG':corr_false_KEGG,
                'REACTOME':corr_false_REACTOME,
                'GWAS':corr_false_GWAS}
            )


    pd.DataFrame(n_FDR).to_csv('random_res/random_results_' + mapmethod + '.csv', index=False)


# Now we can analyse the results


import my_plot_style as mps
f, ax = mps.getplot()

pos = 0
for mapmethod in ('corr', 'TRRUST', 'PPI'):
    res = pd.read_csv('random_res/random_results_' + mapmethod + '.csv')
    print(mapmethod)
    print(res.mean(0))

    ax.bar(range(pos, pos + 4), res.mean(0).values + 1, color=mps.colourscheme())

    pos += 5
ax.set_yscale('log')

labels = []
for _ in range(3):
    for enrtype in res.columns:
        labels.append(enrtype)
    labels.append('')

ax.set_xticklabels(labels, rotation=90)
ax.set_xticks(range(14))
plt.yticks(fontsize=17)
ax.set_ylim([1, 1000])
ax.set_ylabel('annotations log(n + 1)', fontsize=17)
f.savefig('S3.png', bbox_inches='tight')
f.savefig('S3.svg')
