#!/usr/bin/env python
# coding: utf-8

# # load useful libraries

# In[1]:


import scanpy as sc

from warnings import filterwarnings
filterwarnings('ignore')

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# # read data

# In[2]:


adata = sc.read_h5ad('./sc_training.h5ad')
adata


# In[3]:


print(f"There are totally {len(adata.obs.state.unique())} states, they are: {list(adata.obs.state.unique())}.")


# In[4]:


print(f"There are totally {len(adata.obs.condition.unique())} type of gene experiments, which is consistent to description that 66 gene knockouts with unperturb.")


# In[6]:


# transform the data into 
import pandas as pd
df_full = pd.DataFrame(adata.X.toarray())


# ## Challenge 1:
# 
# ### Requirements: 
# 
# Predict the cell-state distribution of the held-out genes: For For each of the 7 held-out knockouts (targeting genes 'Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3'), predict the resulting 5-dimensional vector of cell state proportions (a,b,c,d,e), where ['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other']. 
# > a+b+c+d+e=1

# In[9]:


held_out = ['Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3']
training = list(adata.obs.condition.unique())
for i in held_out:
    if i in training:
        print("This should be in validation dataset.")


# ### Scoring: 
# 
# L1-loss: Mean Absolute Error, example:
# > Actual: (0.2,0.1,0.3,0,0.4)
# 
# > Predict: (0.19,0.11,0.28,0,0.42)
# 
# > Loss = |0.2-0.19|+|0.1-0.11|+|0.3-0.28|+|0-0|+|0.4-0.42|=0.06.

# ## Attempt 1 for Challenge 1: 

# In[ ]:


df_obs = adata.obs.copy()
df_obs = df_obs.reset_index()
df_full.columns = list(adata.var_names)


# ### Inspecting Data: getting some basic information about the data: 

# In[13]:


# 1) why there are only 64 rows in above table? 
experiment_genes = list(df_obs.condition.unique())
expression_genes = list(df_full.columns)
print("Following knockout genes are not features of gene expression matrix:")
for i in experiment_genes:
    if i not in expression_genes and i != "Unperturbed":
        print(i)

print("Following validation knockout genes are not features of gene expression matrix:") 
for i in held_out:
    if i not in expression_genes:
        print(i)


# In[14]:


# change: 
# 2) is there obvious difference between unperturbed gene expression and experimented gene expression? 
# - difference of gene expression for experimenting genes: 
import numpy as np

# getting the average of gene expression: 
unpert_av_expr = np.array(adata.X[adata.obs['condition'] == 'Unperturbed'].mean(axis=0))[0]
df_unpert_av_expr = pd.DataFrame({'gene_name':list(adata.var_names),
                                                  'expr':unpert_av_expr})
knockout_map = {
    'gene_knockout':[],
    'expr':[]
}
for i in range(adata.X.shape[0]):
    knockout = adata.obs.iloc[i]['condition'] # what gene was knocked out 
    if knockout != "Unperturbed":
        df = pd.DataFrame({"gene_name": list(adata.var_names),
                           "expression":adata.X[i].toarray()[0]})
        expr_l = list(df[df['gene_name'] == knockout]['expression'])
        
        if len(expr_l) == 0:
            continue
        knockout_map['gene_knockout'].append(knockout)
        knockout_map['expr'].append(expr_l[0])
        
df_knockout_expr = pd.DataFrame(knockout_map)


# In[15]:


df = df_knockout_expr.merge(df_unpert_av_expr, left_on = "gene_knockout", right_on = "gene_name")
df = df.drop("gene_name", axis = 1)
df.columns = ["gene_knockout", "knockout", "unperturbed"]
df["difference"] = df["knockout"] - df["unperturbed"]
df_sim = df.groupby(["gene_knockout"]).agg({"difference":["mean", "std", "min"],"unperturbed": "mean"})
df_sim.columns = ["mean", "std", "min", "unperturbed"]
df_sim = df_sim.reset_index()


# In[ ]:


df_sim2 = df.merge(df_sim, on = "gene_knockout")
df_sim2["mod_difference"] = df_sim2["difference"] - df_sim2["min"]


# In[18]:


# inspecting the data, the distribution of difference is more likely a gama distribution:
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
def knockout_difference_gamma_distribution(gene, plot = True): 
    dist = df_sim2[df_sim2.gene_knockout == gene]
    try:
        fit_alpha, fit_loc, fit_beta=stats.gamma.fit(dist.difference)
        # plot: 
        if plot: 
            sns.displot(dist, x = "difference", stat="probability")
            start = min(dist.difference)
            end = max(dist.difference)
            x = np.linspace(start, end, 100)
            y = stats.gamma.pdf(x, a = fit_alpha, scale = 1/fit_beta, loc = fit_loc)
            plt.ylim([0,1])
            plt.plot(x,y)
        return {"alpha":fit_alpha, "loc": fit_loc, "scale": 1/fit_beta}
    except:
        pass


# In[19]:


def knockout_difference_beta_distribution(gene, plot = True): 
    dist = df_sim2[df_sim2.gene_knockout == gene]
    try:
        fit_alpha, fit_beta, fit_loc, fit_scale=stats.beta.fit(dist.difference)
        # plot: 
        if plot: 
            sns.displot(dist, x = "difference", stat="probability")
            start = min(dist.difference)
            end = max(dist.difference)
            x = np.linspace(start, end, 100)
            y = stats.beta.pdf(x, a = fit_alpha, b = fit_beta, scale = fit_scale, loc = fit_loc)
            plt.ylim([0,1])
            plt.plot(x,y)
        return {"alpha":fit_alpha, "beta":fit_beta, "loc": fit_loc, "scale": fit_scale}
    except:
        pass


# In[23]:


beta_params = {}

for gene in df_sim.gene_knockout:
    beta_params[gene] = knockout_difference_beta_distribution(gene, plot = False)

beta_params


# In[26]:


import random
def simulate_difference(knockout, samplesize=5000, Knn = 3):
    # sampling cells from dataset:
    df_sample = df_full[controllable_experiment_genes].sample(n = samplesize).transpose()
    distance = []
    # getting 3 nearest neighbours:
    for i in range(len(df_sample)):
        dist = sum(sum(np.square(np.array(df_sample.iloc[i,:] - df_sample[df_sample.index == knockout]))))
        distance.append(dist)

    df_dist = pd.DataFrame({"gene": controllable_experiment_genes,
                            "distance": distance})

    nn_genes = df_dist.sort_values(by = ["distance"])[1:Knn+1].reset_index()
    total_dist = sum(nn_genes.distance)
    total_diff = 0
    unpert_knock_ave = np.array(df_sim.unperturbed[df_sim.gene_knockout == knockout])
    
    for i in range(Knn):
        # total_diff += sum(df_comparison[df_comparison.index == nn_genes.gene[i]].difference*nn_genes.distance[i])
        params = beta_params[nn_genes.gene[i]]
        cur = nn_genes.gene[i]
        unpert_cur_ave = np.array(df_sim.unperturbed[df_sim.gene_knockout == cur])
        
        a = params["alpha"]
        b = params["beta"]
        loc = params["loc"]
        scale= params["scale"]
        s = stats.beta.rvs(a, b, loc, scale, size=1)
        
        total_diff += (s*(unpert_knock_ave+1)/(unpert_cur_ave+1))*nn_genes.distance[i]
        
    sim = sum(total_diff/total_dist)
    
    return sim, nn_genes


# In[114]:


experiment_genes = list(df_obs.condition.unique())
experiment_genes = list(set(experiment_genes).difference({"Unperturbed", "Fzd1", "P2rx7"}))
experiment_genes.extend(['Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3'])
len(experiment_genes)

validation_test = ['Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3']
validation = ['Aqr', 'Bach2', 'Bhlhe40']
test = ['Ets1', 'Fosb', 'Mafk', 'Stat3']
controllable_experiment_genes = experiment_genes.copy()


# In[111]:


def simulate_state(knockout, samplesize=5000, Knn = 3, echos = 10):
    # sampling cells from dataset:
    df_sample = df_full[controllable_experiment_genes].sample(n = samplesize).transpose()
    distance = []
    # getting 3 nearest neighbours:
    for i in range(len(df_sample)):
        dist = sum(sum(np.square(np.array(df_sample.iloc[i,:] - df_sample[df_sample.index == knockout]))))
        distance.append(dist)

    df_dist = pd.DataFrame({"gene": controllable_experiment_genes,
                            "distance": distance})

    nn_genes = df_dist.sort_values(by = ["distance"]).reset_index()
    
    total_weight = 0 
    state_ratio = 0
    var = nn_genes.distance.var()+10
    
    for e in range(echos): 
        k = 0
        i = 0 
    
        while k < Knn and i < len(nn_genes):
            cur = nn_genes.gene[i]
            if cur not in validation_test and cur != knockout:
            
                k += 1
        
                weight = np.exp(-nn_genes.distance[i]/(2*var))
                total_weight += weight
            
                state_count = np.array(df_obs[df_obs.condition == cur].groupby("state").state.count())
                state_ratio += state_count*weight/sum(state_count)
            i = i+1
        
        
    
    state_ratio = state_ratio/total_weight
    state_name = ["cycling", "effector", "other", "progenitor", "terminal exhausted"]
    state_map = {}
    
    for i in range(5):
        state_map[state_name[i]] = state_ratio[i]
        
    return state_map


# In[122]:


validation_map = {}
for v in validation:
    vmap = simulate_state(knockout = v, samplesize=5000, Knn = 3)
    formal_order = ['progenitor','effector','terminal exhausted','cycling','other']
    vlist = []
    for i in range(5):
        vlist.append(vmap[formal_order[i]])
    validation_map[v] = vlist 
        
validation_output = pd.DataFrame(validation_map).transpose()
validation_output.columns = formal_order
validation_output


# In[125]:


#validation_output.to_csv("validation_output.csv")


# In[121]:


test_map = {}
for t in test:
    tmap = simulate_state(knockout = t, samplesize=5000, Knn = 3)
    formal_order = ['progenitor','effector','terminal exhausted','cycling','other']
    tlist = []
    for i in range(5):
        tlist.append(tmap[formal_order[i]])
    test_map[t] = tlist 
        
test_output = pd.DataFrame(test_map).transpose()
test_output.columns = formal_order
test_output


# In[124]:


#test_output.to_csv("test_output.csv")

