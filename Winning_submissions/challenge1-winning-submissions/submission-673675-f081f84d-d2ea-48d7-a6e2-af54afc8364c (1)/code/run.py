import os

os.system('wget https://saturn-public-data.s3.us-east-2.amazonaws.com/cancer-immunotherapy-challenge/data/sc_training.h5ad')
os.system('pip install scanpy')

import scanpy 
import pandas as pd
import numpy as np
import scanpy 
from sklearn.neighbors import NearestNeighbors


adata = scanpy.read_h5ad('sc_training.h5ad')
counts = adata.obs.condition.value_counts()
good_genes = counts[counts > 30].index

nn_genes = set(good_genes) - set(['Fzd1', 'P2rx7', 'Unperturbed'])
nn_genes = list(nn_genes) + [ 'Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3']

mean_data = {}
targets = ['progenitor', 'effector',  'terminal exhausted' ,'cycling', 'other']
for i, gene in enumerate(nn_genes[:-7]):
    for target in targets:
        sh = adata[adata.obs.condition == gene].shape
        mean_ = 0
        if len(sh):
            sh1 = adata[(adata.obs.condition == gene) & (adata.obs.state == target)].shape
            if len(sh1):
                mean_ = sh1[0] / sh[0]
        mean_data[(i, target)] = mean_

clf = NearestNeighbors(n_neighbors = 5, p = 1)

keys = []
for gene in nn_genes:
    keys += [np.where(adata.var.index == gene)[0][0]]

clf.fit(adata.X.T[keys[:-7]])
answer = []
for ind in clf.kneighbors(adata.X.T[keys[-7:]])[1]:
    tmp_answer = []
    for target in targets:
        tmp_answer += [np.mean([mean_data[(i, target)] for i in ind])]
    answer += [tmp_answer]

predict_df = pd.DataFrame(answer[:3], columns = ['a_i','b_i','c_i','d_i','e_i'])
predict_df['gene'] = ['Aqr', 'Bach2', 'Bhlhe40']
predict_df = predict_df[['gene','a_i','b_i','c_i','d_i','e_i']]
predict_df.to_csv('validation_output.csv', index = None)

predict_df = pd.DataFrame(answer[3:], columns = ['a_i','b_i','c_i','d_i','e_i'])
predict_df['gene'] = ['Ets1', 'Fosb', 'Mafk', 'Stat3']
predict_df = predict_df[['gene','a_i','b_i','c_i','d_i','e_i']]
predict_df.to_csv('test_output.csv', index = None)