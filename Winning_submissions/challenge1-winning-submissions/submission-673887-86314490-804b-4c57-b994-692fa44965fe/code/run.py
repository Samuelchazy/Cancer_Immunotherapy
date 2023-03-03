def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

seed_value = 34
import random
random.seed(seed_value)
import numpy as np
np.random.seed(seed_value)
import pandas as pd
from sklearn.neighbors import KNeighborsRegressor

import scanpy as sc
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor='red')

adata = sc.read_h5ad('input/sc_training.h5ad')
adata

adata.obs

adata.obs["RNA"] = adata.obs.index
adata_arr = adata.obs.values

adata_arr

def extract_rna(brc):
    dx = 0
    while dx<len(brc) and brc[dx]!="_":
        dx += 1
    dx+=1
    rna=""
    while dx<len(brc) and brc[dx]!="-":
        rna = rna + str(brc[dx])
        dx += 1
    return rna

condtions = []
barcode_states = {}
states_map = {}
rev_states_map = {}
for rw in adata_arr:
    barcode_states[extract_rna(rw[4])] = rw[1]
    if rw[1] not in states_map:
        rev_states_map[len(states_map)] = rw[1]
        states_map[rw[1]] = len(states_map)
    if rw[2] not in condtions:
        condtions.append(rw[2])

targets1 = ["Aqr","Bach2","Bhlhe40","Ets1","Fosb","Mafk","Stat3"]

all_grna = {}
for el in targets1:
    all_grna[el] = len(all_grna)
for el in condtions:
    if el not in all_grna:
        all_grna[el] = len(all_grna)

clon_info = pd.read_csv('input/clone_information.csv')
clon_info

clon_info=clon_info.values
clon_info

numbers = {}

for rw in clon_info:
    wh = ""
    for key in all_grna:
        if key in rw[1]:
            if len(wh)<len(key):
                wh = key
    if wh!="":
        if wh not in numbers:
            numbers[wh] = 0.0
        numbers[wh] += 1.0

types = {}
for kn in condtions:
    types[kn] = {}
    for st in states_map:
        types[kn][st] = 0.0

for i in range(0,len(adata_arr)):
    kn = adata_arr[i][2]
    st = adata_arr[i][1]
    types[kn][st] += 1.0

for kn in types:
    summ = 0.0
    for st in types[kn]:
        summ += types[kn][st]
    for st in types[kn]:
        types[kn][st] /= summ

states = ["progenitor","effector","terminal exhausted","cycling","other"]

statesx = []
progenitory, effectory, terminaly, cyclingy, othery = [], [], [], [], []
states_test = []

for cnd in targets1:
    states_test.append(numbers[cnd])

for cnd in condtions:
    if (cnd in numbers) and (cnd in types) and (cnd not in targets1):
        statesx.append(numbers[cnd])
        progenitory.append(types[cnd]["progenitor"])
        effectory.append(types[cnd]["effector"])
        terminaly.append(types[cnd]["terminal exhausted"])
        cyclingy.append(types[cnd]["cycling"])
        othery.append(types[cnd]["other"])

statesx = pd.DataFrame(statesx)
states_test = pd.DataFrame(states_test)
progenitory = pd.DataFrame(progenitory)
effectory = pd.DataFrame(effectory)
terminaly = pd.DataFrame(terminaly)
cyclingy = pd.DataFrame(cyclingy)
othery = pd.DataFrame(othery)

predicts = []

knn = KNeighborsRegressor(n_neighbors=15,weights='uniform')
knn.fit(statesx,progenitory)

predicts.append(knn.predict(states_test))

knn = KNeighborsRegressor(n_neighbors=15,weights='uniform')
knn.fit(statesx,effectory)

predicts.append(knn.predict(states_test))

knn = KNeighborsRegressor(n_neighbors=15,weights='uniform')
knn.fit(statesx,terminaly)

predicts.append(knn.predict(states_test))

knn = KNeighborsRegressor(n_neighbors=5,weights='uniform')
knn.fit(statesx,cyclingy)

predicts.append(knn.predict(states_test))

knn = KNeighborsRegressor(n_neighbors=15,weights='uniform')
knn.fit(statesx,othery)

predicts.append(knn.predict(states_test))

pred = []
for i in range(len(predicts[0])):
    rw = []
    for j in range(len(predicts)):
        rw.append(predicts[j][i][0])
    summ = sum(rw)
    for j in range(0,len(rw)):
        rw[j] = rw[j] / summ
    pred.append(rw)

knock = {}
for kn in targets1:
    knock[kn] = {}
    for st in states_map:
        knock[kn][st] = 0.0

for i in range(0,len(pred)):
    for j in range(5):
        knock[targets1[i]][states[j]] += pred[i][j]

for kn in targets1:
    summ = 0.0
    for st in states_map:
        summ += knock[kn][st]
    for st in states_map:
        knock[kn][st] /= summ

validation = []
for kn in ["Aqr","Bach2","Bhlhe40"]:
    rw = [kn]
    for st in ["progenitor","effector","terminal exhausted","cycling","other"]:
        rw.append(knock[kn][st])
    validation.append(rw)

test = []
for kn in ["Ets1","Fosb","Mafk","Stat3"]:
    rw = [kn]
    for st in ["progenitor","effector","terminal exhausted","cycling","other"]:
        rw.append(knock[kn][st])
    test.append(rw)

validation = pd.DataFrame(validation,columns=["gene","a_i","b_i","c_i","d_i","e_i"])
test = pd.DataFrame(test,columns=["gene","a_i","b_i","c_i","d_i","e_i"])

validation.to_csv("../solution/validation_output.csv",index=False)
test.to_csv("../solution/test_output.csv",index=False)