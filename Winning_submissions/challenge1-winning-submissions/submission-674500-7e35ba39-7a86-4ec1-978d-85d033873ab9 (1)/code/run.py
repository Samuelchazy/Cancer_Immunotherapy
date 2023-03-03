# %%
# Import essential libraries
import scanpy as sc
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from tqdm.auto import tqdm 

import pickle
from statistics import mean

# Import PyTorch
import torch
from torch import nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

# Import Fast.AI tabular
from fastai.tabular.all import *
from fastai.imports import *
from fastai.torch_core import *
from fastai.learner import *

# Set up device agnostic code
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f'Device used: {device}')

from warnings import filterwarnings
filterwarnings('ignore')

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# %%
adata = sc.read_h5ad('/home/kokyriakidis/Downloads/topcoder/sc_training.h5ad')

with open('/home/kokyriakidis/Downloads/states_dictionary.pkl', 'rb') as f:
    states_dictionary = pickle.load(f)

#reactome_table = pd.read_csv('/home/kokyriakidis/Downloads/Reactome_2022_table_all.txt', delimiter='\t')

with open('/home/kokyriakidis/Downloads/diff_expression_double_shrunk_lfc2_001.pkl', 'rb') as handle:
    diff_expression = pickle.load(handle)

# Get cell names
cell_names = adata.obs_names

# Get gene names|
gene_names = adata.var_names

len(cell_names), len(gene_names)

# Validation and Test genes
val_test_genes = ['Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3']

# Create metadata dataframe
metadata = adata.obs[['state','condition']]
#metadata = pd.read_csv('metadata.csv', index_col = 0)

# Initiate sparse matrix
sparse_matrix = adata.X

# (28697, 15077) -> (cells, genes) (rows, cols)
sparse_matrix_shape = sparse_matrix.toarray().shape
print(f'There are {sparse_matrix_shape[0]} cells (rows).')
print(f'There are {sparse_matrix_shape[1]} genes (columns).')

# Sparse normalized dataframe
sparse_df = pd.DataFrame(sparse_matrix.toarray())
# Sparse row dataframe
#sparse_df = pd.DataFrame(adata.layers['rawcounts'].toarray())
sparse_df.index = cell_names
sparse_df.columns = gene_names

print(f'There are {sparse_df.shape[0]} cells (rows).')
print(f'There are {sparse_df.shape[1]} genes (columns).')

with open('/home/kokyriakidis/Downloads/diff_expression_double_shrunk_lfc2_001.pkl', 'rb') as handle:
    diff_expression = pickle.load(handle)

# Get cell names
cell_names = adata.obs_names

# Get gene names|
gene_names = adata.var_names

len(cell_names), len(gene_names)

# Validation and Test genes
val_test_genes = ['Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3']

# Create metadata dataframe
metadata = adata.obs[['state','condition']]
#metadata = pd.read_csv('metadata.csv', index_col = 0)

# Initiate sparse matrix
sparse_matrix = adata.X

# (28697, 15077) -> (cells, genes) (rows, cols)
sparse_matrix_shape = sparse_matrix.toarray().shape
print(f'There are {sparse_matrix_shape[0]} cells (rows).')
print(f'There are {sparse_matrix_shape[1]} genes (columns).')

# Sparse normalized dataframe
sparse_df = pd.DataFrame(sparse_matrix.toarray())
# Sparse row dataframe
#sparse_df = pd.DataFrame(adata.layers['rawcounts'].toarray())
sparse_df.index = cell_names
sparse_df.columns = gene_names

print(f'There are {sparse_df.shape[0]} cells (rows).')
print(f'There are {sparse_df.shape[1]} genes (columns).')

# %%
# Find the T-cell distributions for each knocked out gene
distributions_list = []
conditions = list(metadata.condition)
distributions_dict = {}
sum_dict = {}
sum_dict_2 = {}
for con in set(conditions):
    sub_metadata = metadata[metadata['condition'] == con]
    counts = sub_metadata['state'].value_counts()
    sum = counts.sum()
    sum_dict_2[con] = sum
    distribution = counts/sum
    try:
        cycling = distribution['cycling']
        sum_cycling = counts['cycling']
    except:
        cycling = 0
        sum_cycling=0
        pass

    try:
        t_exh = distribution['terminal exhausted']
        sum_t_exh = counts['terminal exhausted']
    except:
        t_exh = 0
        sum_t_exh=0
        pass

    try:
        effector = distribution['effector']
        sum_effector = counts['effector']
    except:
        effector = 0
        sum_effector=0
        pass
    
    try:
        progenitor = distribution['progenitor']
        sum_progenitor = counts['progenitor']
    except:
        progenitor = 0
        sum_progenitor=0
        pass

    try:
        other = distribution['other']
        sum_other = counts['other']
    except:
        other = 0
        sum_other=0
        pass

    row = [progenitor, effector, t_exh, cycling, other]
    sum_row = [sum_progenitor, sum_effector, sum_t_exh, sum_cycling, sum_other]
    distributions_dict[con] = row
    sum_dict[con] = sum_row
for con in conditions:
    row = distributions_dict[con]
    distributions_list.append(row)
distributions_df = pd.DataFrame(distributions_list, columns = ['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other'])
distributions_df.index = cell_names

# Add the distribution columns to the metadata dataframe
metadata = pd.concat([metadata, distributions_df], axis=1)


# %%
delete_genes = []
for con in sum_dict_2.keys():
    if distributions_dict[con] == [1,0,0,0,0]:
        #print(f'Gene: {con}\t| Cell count: {sum_dict_2[con]}-{sum_dict[con]} \t| Distribution: {distributions_dict[con]}')
        delete_genes.append(con)
    elif sum_dict_2[con] < 11: 
        #print(f'Gene: {con}\t| Cell count: {sum_dict_2[con]}-{sum_dict[con]} \t| Distribution: {distributions_dict[con]}')
        delete_genes.append(con)



# %%
delete_cells_metadata = list(metadata[metadata['condition'].isin(delete_genes)].index)
print(len(set(metadata.condition)))
metadata = metadata.drop(delete_cells_metadata)
sparse_df = sparse_df.drop(delete_cells_metadata)
#sparse_df = sparse_df.drop(delete_genes, axis=1)

#print(f'There are {sparse_df.shape[0]} cells (rows).')
#print(f'There are {sparse_df.shape[1]} genes (columns).')

# %%
all_genes = []
# iterate over files in
# that directory
for key in  diff_expression.keys():
    stats_df = diff_expression[key]
    all_genes = list(all_genes) + list(stats_df.index)

cols = list(set(all_genes))
len(all_genes), len(set(all_genes))

# %%
states = ['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other']
state = 'progenitor'
indexes = {}
for state in states:
    i = list(metadata[(metadata['state']==state) & (metadata[state] == metadata[states].max(axis=1))].index) #
    indexes[state] = i
    print(f'State {state} has {len(i)} rows')

# %% [markdown]
# ### Set up training and validation datasets

# %%
df = sparse_df[cols]
df_nozeros = df.copy()

# %%
# replace zero cell values with mean values
states = ['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other']

for state in tqdm(states):
    # Find indexes of cells in state
    cells_index = metadata[(metadata['state']==state)].index
    state_dict = states_dictionary[state]
    sub_df = df_nozeros.loc[cells_index]
    for gene in tqdm(state_dict.keys(), leave=False):
        if gene in df_nozeros.columns:
            i = sub_df[sub_df[gene] == 0].index
            df_nozeros.loc[i, gene] = state_dict[gene]

# %%
len(cols)

# %%
knockout_genes = list(set(metadata.condition))
knockout_genes.remove('Unperturbed')
knockout_genes.remove('P2rx7')
knockout_genes.remove('Fzd1')

# %%
genes_of_interest = knockout_genes + val_test_genes

# %%
# # How many paths to keep
# NUM_PATHS = 3
# paths_genes = {}
# for gene in genes_of_interest:
#     sub_df = reactome_table[reactome_table['Genes'].str.contains(gene.upper())]
#     sub_df = sub_df[sub_df['Genes'].eq(gene.upper())==False].reset_index(drop=True)
#     if sub_df.empty:
#         continue
#     else:
#         temp_1 = [gene] + val_test_genes
#         temp_1 = list(map(str.lower,temp_1))
#         delete_index = []
#         for index, row in sub_df.iterrows():
#             temp_2 = row['Genes'].split(';')
#             temp_2 = list(map(str.lower,temp_2))
#             if (all(x in temp_1 for x in temp_2)):
#                 delete_index.append(index)
#         sub_df = sub_df.drop(delete_index)
#         if sub_df.empty:
#             continue
#         else:
#             sub_df = sub_df.nlargest(NUM_PATHS,['Combined Score'])
#             genes_str = ';'.join(list(sub_df['Genes']))
#             genes_list = list(set(genes_str.split(';')))
#             for val_test_gene in val_test_genes:
#                 if val_test_gene in genes_list:
#                     genes_list.remove(val_test_gene)
#             genes_list.remove(gene.upper())
#             paths_genes[gene] = [item.capitalize() for item in genes_list]

# %%
# Remove knockout genes
for gene in knockout_genes:
    if gene in cols:
        df_nozeros = df_nozeros.drop(gene, axis=1)
        print(f'{gene} dropped')

# %%
# create validation (save_df) and training set (df_filtered)
validation_genes = [ 'Tox2', 'Tcf3', 'Hmgb2', 'Dvl3', 'Tpt1', 'Ctnnb1','Foxm1']

#validation_genes =  ['Tox2', 'Irf2', 'Zeb2']
#validation_genes =  ['Foxm1', 'Myb', 'Ctnnb1']
true_values = pd.DataFrame(columns=['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other'])
delete_cells = []

for gene in validation_genes:
    sub_df = metadata[metadata['condition'] == gene]
    distr = [sub_df['progenitor'][0], sub_df['effector'][0], sub_df['terminal exhausted'][0], sub_df['cycling'][0], sub_df['other'][0]]
    true_values.loc[gene]=distr
    index = sub_df.index
    delete_cells += list(index)
save_df = df_nozeros.loc[delete_cells]
df_filtered = df_nozeros.drop(delete_cells)
print(df_filtered.shape, save_df.shape)
#print(f'There are {sparse_df.shape[0]} cells (rows).')
#print(f'There are {sparse_df.shape[1]} genes (columns).')

display(true_values)

# %%
# # Add weighted mean lines to the training set
# distr = {}
# for gene in knockout_genes:
#     if gene in validation_genes: # Skip if the gene belongs to the extracted validation set
#         continue
#     p = metadata[(metadata['condition']==gene)]['progenitor'].mean()
#     e = metadata[(metadata['condition']==gene)]['effector'].mean()
#     t = metadata[(metadata['condition']==gene)]['terminal exhausted'].mean()
#     c = metadata[(metadata['condition']==gene)]['cycling'].mean()
#     distr[gene] = [p, e, t, c]


# for i, gene in enumerate(distr.keys()):
#     d = distr[gene]
#     cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='progenitor')].index
#     #sub_p = df_filtered.loc[cell_index].mean() * d[0]
#     sub_p = df_filtered.loc[cell_index].median()

#     cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='effector')].index
#     #sub_e = df_filtered.loc[cell_index].mean()* d[1]
#     sub_e = df_filtered.loc[cell_index].median()

#     cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='terminal exhausted')].index
#     #sub_t = df_filtered.loc[cell_index].mean()* d[2]
#     sub_t = df_filtered.loc[cell_index].median()

#     cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='cycling')].index
#     #sub_c = df_filtered.loc[cell_index].mean()* d[3]
#     sub_c = df_filtered.loc[cell_index].median()

#     f = (sub_p + sub_e + sub_t + sub_c) / (np.sum(d[0:-1]))
#     df_filtered.loc[f'{gene}_weighted_mean'] = f



# %%
cols = df_filtered.columns
df_mean = pd.DataFrame(columns=cols)
states = ['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other']
for gene in set(metadata.condition):
   if gene in validation_genes:
       continue
   rows = metadata[(metadata['condition']==gene)].index
   row = df_filtered.loc[rows]
   df_mean.loc[f'{gene}_mean'] = row.mean()
df_mean

# %%
# Reduced rows
index_1 = list(set(indexes['progenitor']) - set(delete_cells)) 
index_2 = random.sample(list(set(indexes['effector']) - set(delete_cells)), k=300)
index_3 = list(set(indexes['terminal exhausted']) - set(delete_cells)) #random.choices(list(set(indexes['terminal exhausted']) - set(delete_cells)), k=400) #list(set(indexes['terminal exhausted']) - set(delete_cells))
index_4 = random.sample(list(set(indexes['cycling']) - set(delete_cells)), k=300)
index_5 = list(set(indexes['other']) - set(delete_cells))

df_f = df_filtered.loc[index_1 + index_2 + index_3 +  index_4 + index_5]

# Add target
df_f['progenitor'] = metadata.loc[df_f.index]["progenitor"]
df_f['effector'] = metadata.loc[df_f.index]["effector"]
df_f['terminal exhausted'] = metadata.loc[df_f.index]["terminal exhausted"]
df_f['cycling'] = metadata.loc[df_f.index]["cycling"]
df_f['other'] = metadata.loc[df_f.index]["other"]

# %%
len(index_2)

# %% [markdown]
# ### Model training

# %%
#torch.manual_seed(0)
dls = TabularDataLoaders.from_df(df_f, 
                                 #procs=procs,
                                 y_names=["progenitor","effector","terminal exhausted", "cycling", "other"], 
                                 shuffle=True,
                                 val_shuffle=True,
                                 device=device)
#torch.manual_seed(0)
# Initiate tabular learner
metrics=[mse, rmse, mae, R2Score()]
#torch.manual_seed(0)
#learn = tabular_learner(dls,  metrics=metrics, y_range = (0,1), loss_func=L1LossFlat(reduction='sum'))

learn = tabular_learner(dls, metrics=metrics, y_range = (0,1), loss_func=nn.MSELoss(reduction='sum'))

# Find optimum learning rate
#torch.manual_seed(0)
lrs = learn.lr_find(suggest_funcs=(minimum, steep, valley, slide))

# Train model
#torch.manual_seed(0)
test = learn.fit_one_cycle(50, lrs.valley, cbs=[EarlyStoppingCallback(monitor='mae', comp=np.less, min_delta=0.0001, patience=50)])

# %% [markdown]
# ### Model evaluation

# %%
def calculate_mae(act, pred):
    sums = 0
    for index in act.index:
        y_true = act.loc[index].values
        predictions = pred.loc[index].values
        y_true, predictions = np.array(y_true), np.array(predictions)
        sum = np.sum(np.abs(y_true - predictions))
        sums += sum

        print(sum)
    print(f'MAE: {sums/len(act.index)}')
    #return sums/7

# %%
distr = {}
for gene in validation_genes:
    p = metadata[(metadata['condition']==gene)]['progenitor'].mean()
    e = metadata[(metadata['condition']==gene)]['effector'].mean()
    t = metadata[(metadata['condition']==gene)]['terminal exhausted'].mean()
    c = metadata[(metadata['condition']==gene)]['cycling'].mean()
    distr[gene] = [p, e, t, c]

# %%
final_test = pd.DataFrame(columns=save_df.columns)
for i, gene in enumerate(distr.keys()):
    d = distr[gene]
    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='progenitor')].index
    sub_p = save_df.loc[cell_index].mean() * d[0]

    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='effector')].index
    sub_e = save_df.loc[cell_index].mean()* d[1]

    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='terminal exhausted')].index
    sub_t = save_df.loc[cell_index].mean()* d[2]

    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='cycling')].index
    sub_c = save_df.loc[cell_index].mean()* d[3]

    f = (sub_p + sub_e + sub_t + sub_c) / (np.sum(d[0:-1]))
    final_test.loc[gene] = f

final_test

# %%
test_dl = learn.dls.test_dl(final_test, num_workers=0, shuffle=False)
preds = learn.get_preds(dl=test_dl)

clas_1 = preds[0][0]
clas_2 = preds[0][1]
clas_3 = preds[0][2]
clas_4 = preds[0][3]
clas_5 = preds[0][4]
clas_6 = preds[0][5]
clas_7 = preds[0][6]

validation_list = [clas_1.tolist(), clas_2.tolist(), clas_3.tolist(), clas_4.tolist(), clas_5.tolist(), clas_6.tolist(), clas_7.tolist()]
validation_set = pd.DataFrame(validation_list, columns=['a_i','b_i','c_i','d_i','e_i'])
validation_set.insert(0, 'genes', final_test.index)
validation_set = validation_set.set_index('genes')
sums = validation_set.sum(axis=1)
for index in sums.index:
    sum = sums[index]
    validation_set.loc[index] = validation_set.loc[index]/sum
validation_set

# %%
true_values

# %%
calculate_mae(validation_set, true_values)

# %% [markdown]
# ## Use best model

# %%
#learn.export('/home/kokyriakidis/Downloads/mean_mae_01909_no_knockouts.pkl')
learn = load_learner("/home/kokyriakidis/Downloads/mean_mae_01909_no_knockouts.pkl")

# %%
# Testing for Hmgb2
gene_lists_1 = 'HMGB1' #
#gene_lists_1 = 'SATB1;HMGB1;CTNNB1'
#gene_lists_1 = paths_genes["HMGB2"]
#gene_lists_1 = 'HMGB1;CTNNB1;IRF2;RELA;ATF2;NR3C1;LEF1;EZH2'
# if ";" in gene_lists_1:
#     gene_list_2 = gene_lists_1.split(';')
# else:
#     gene_list_2 = paths_genes["HMGB1"]
try:
    gene_list_2 = paths_genes["HMGB1"]
except:
    try:
        gene_list_2 = gene_lists_1.split(';')
    except:
        gene_list_2 = gene_lists_1
gene_list=[]
for g in gene_list_2:
    gene_list.append(g.capitalize())
distr = {}
for gene in gene_list:
    if gene not in knockout_genes:
        continue
    if gene in val_test_genes or gene == 'Hmgb2':
        continue
    p = metadata[(metadata['condition']==gene)]['progenitor'].mean()
    e = metadata[(metadata['condition']==gene)]['effector'].mean()
    t = metadata[(metadata['condition']==gene)]['terminal exhausted'].mean()
    c = metadata[(metadata['condition']==gene)]['cycling'].mean()
    distr[gene] = [p, e, t, c]

# %%
final_df = pd.DataFrame(columns=df_nozeros.columns)
for i, gene in enumerate(distr.keys()):
    d = distr[gene]
    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='progenitor')].index
    sub_p = df_nozeros.loc[cell_index].mean() * d[0]

    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='effector')].index
    sub_e = df_nozeros.loc[cell_index].mean()* d[1]

    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='terminal exhausted')].index
    sub_t = df_nozeros.loc[cell_index].mean()* d[2]

    cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='cycling')].index
    sub_c = df_nozeros.loc[cell_index].mean()* d[3]

    f = (sub_p + sub_e + sub_t + sub_c) / (np.sum(d[0:-1]))
    final_df.loc[gene] = f

validation_rows = pd.DataFrame(columns=df_nozeros.columns)
validation_rows.loc['Hmgb2'] = final_df.mean()



# %%
test_dl = learn.dls.test_dl(validation_rows, num_workers=0, shuffle=False)
preds = learn.get_preds(dl=test_dl)

clas_1 = preds[0][0]

validation_list = [clas_1.tolist()]
validation_set = pd.DataFrame(validation_list, columns=['a_i','b_i','c_i','d_i','e_i'])
validation_set.insert(0, 'genes', validation_rows.index)
validation_set = validation_set.set_index('genes')
sums = validation_set.sum(axis=1)
for index in sums.index:
    sum = sums[index]
    validation_set.loc[index] = validation_set.loc[index]/sum
validation_set

# %%
true_values

# %%
calculate_mae(validation_set, true_values)

# %% [markdown]
# ### Read STRING database

# %%
variant_file = open("/home/kokyriakidis/Downloads/10090.protein.info.v11.5.txt")
read_variant_file = csv.reader(variant_file, delimiter="\t")

vfile = []

for row in read_variant_file:
    vfile.append(row)

header = vfile[0]

body = vfile[1:]

ensmus_to_gene_name = {}

for row in range(0, len(body)):
    if body[row][0] not in ensmus_to_gene_name:
        ensmus_to_gene_name[body[row][0]] = body[row][1]

ensmus_to_gene_name


# %%
variant_file = open("/home/kokyriakidis/Downloads/10090.protein.links.v11.5.txt")
read_variant_file = csv.reader(variant_file, delimiter=" ")

vfile = []

for row in read_variant_file:
    vfile.append(row)

header = vfile[0]

body = vfile[1:]

string_mus_db = []

body[0]

for row in range(0, len(body)):
    if (ensmus_to_gene_name[body[row][0]] in genes_of_interest) and (ensmus_to_gene_name[body[row][1]] in genes_of_interest):
        string_mus_db.append([ensmus_to_gene_name[body[row][0]], ensmus_to_gene_name[body[row][1]], int(body[row][2])])

string_mus_db = sorted(string_mus_db, key=lambda x: (x[2], x[0]), reverse=True)

#string_mus_db

# %%
string_mus_db_associations_of_interest = {}

for association in range(0, len(string_mus_db)):
    if string_mus_db[association][0] not in string_mus_db_associations_of_interest:
        string_mus_db_associations_of_interest[string_mus_db[association][0]] = [[], []]
        string_mus_db_associations_of_interest[string_mus_db[association][0]][0].append(string_mus_db[association][1])
        string_mus_db_associations_of_interest[string_mus_db[association][0]][1].append(string_mus_db[association][2])
    if (string_mus_db[association][2] >= 600) and (string_mus_db[association][1] not in string_mus_db_associations_of_interest[string_mus_db[association][0]][0]):
        string_mus_db_associations_of_interest[string_mus_db[association][0]][0].append(string_mus_db[association][1])
        string_mus_db_associations_of_interest[string_mus_db[association][0]][1].append(string_mus_db[association][2])

#string_mus_db_associations_of_interest

# %%
string_mus_db_associations_of_interest_final = {}

for key in string_mus_db_associations_of_interest:
    if len(string_mus_db_associations_of_interest[key][0]) == 1 and (len(set(string_mus_db_associations_of_interest[key][0]).intersection(val_test_genes)) == len(string_mus_db_associations_of_interest[key][0])):
        string_mus_db_associations_of_interest_final[key] = string_mus_db_associations_of_interest[string_mus_db_associations_of_interest[key][0][0]][0]
    if (len(string_mus_db_associations_of_interest[key][0]) < 4) and (len(string_mus_db_associations_of_interest[key][0]) > 1) and (len(set(string_mus_db_associations_of_interest[key][0]).intersection(val_test_genes)) == len(string_mus_db_associations_of_interest[key][0])):
        string_mus_db_associations_of_interest_final[key] = string_mus_db_associations_of_interest[string_mus_db_associations_of_interest[key][0][0]][0]
    if len(string_mus_db_associations_of_interest[key][0]) < 4 and (len(set(string_mus_db_associations_of_interest[key][0]).intersection(val_test_genes)) < len(string_mus_db_associations_of_interest[key][0])):
        string_mus_db_associations_of_interest_final[key] = string_mus_db_associations_of_interest[key][0]
    if (len(string_mus_db_associations_of_interest[key][0]) >= 4) and (len(set(string_mus_db_associations_of_interest[key][0]).intersection(val_test_genes)) < len(string_mus_db_associations_of_interest[key][0])):
        index = 0
        for i in range(3, len(string_mus_db_associations_of_interest[key][0])):
            if string_mus_db_associations_of_interest[key][1][i] < 950:
                index = i
                break
        string_mus_db_associations_of_interest_final[key] = string_mus_db_associations_of_interest[key][0][0:index]
    


# %%
string_mus_db_associations_of_interest_final["Aqr"], string_mus_db_associations_of_interest_final["Bach2"], string_mus_db_associations_of_interest_final["Bhlhe40"], string_mus_db_associations_of_interest_final["Ets1"], string_mus_db_associations_of_interest_final["Fosb"], string_mus_db_associations_of_interest_final["Mafk"], string_mus_db_associations_of_interest_final["Stat3"]

# %% [markdown]
# ## Generate validation and test files

# %%
# relative_genes = {'Aqr': 'YY1;RELA',
#                   'Bach2': 'LRP1',
#                   'Bhlhe40': 'CREM;NR3C1;HIF1A',
#                   'Ets1': 'ATF2;SATB1;LEF1;TCF7;HMGB2;CREM;PRDM1;NR3C1;ETS1;IKZF3;LITAF;HIF1A;BACH2;RELA;FOXO1;YY1;MYB;STAT4;SOX4;EOMES;EGR1;STAT3;RUNX2;FOXP1;NR4A2;NR4A1;ZEB2;ELF1;NR4A3;IRF2;BHLHE40;FOSB;MAFK;TCF3;IRF9',
#                   'Fosb': 'YY1',
#                   'Mafk':'MYB;IRF2',
#                   'Stat3':'STAT4;IL12RB1'}


relative_genes = {'Aqr': string_mus_db_associations_of_interest_final['Aqr'],
                  'Bach2': string_mus_db_associations_of_interest_final['Bach2'],
                  'Bhlhe40': string_mus_db_associations_of_interest_final['Bhlhe40'],
                  'Ets1': string_mus_db_associations_of_interest_final['Ets1'],
                  'Fosb': string_mus_db_associations_of_interest_final['Fosb'],
                  'Mafk': string_mus_db_associations_of_interest_final['Mafk'],
                  'Stat3': string_mus_db_associations_of_interest_final['Stat3'] }

validation_test_df = pd.DataFrame(columns=df_nozeros.columns)

for key in relative_genes.keys():
    temp = pd.DataFrame(columns=df_nozeros.columns)

    distr={}
    gene_lists_1 = relative_genes[key]
    try:
        gene_list_2 = gene_lists_1.split(';')
    except:
        gene_list_2 = gene_lists_1
    gene_list=[]
    for g in gene_list_2:
        gene_list.append(g.capitalize())
    #print(gene_list)
    for gene in gene_list:
        if gene not in knockout_genes:
            continue
        if gene in val_test_genes or gene == key:
            continue
        p = metadata[(metadata['condition']==gene)]['progenitor'].mean()
        e = metadata[(metadata['condition']==gene)]['effector'].mean()
        t = metadata[(metadata['condition']==gene)]['terminal exhausted'].mean()
        c = metadata[(metadata['condition']==gene)]['cycling'].mean()
        distr[gene] = [p, e, t, c]
    
    for i, gene in enumerate(distr.keys()):
        d = distr[gene]
        cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='progenitor')].index
        sub_p = df_nozeros.loc[cell_index].mean() * d[0]

        cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='effector')].index
        sub_e = df_nozeros.loc[cell_index].mean()* d[1]

        cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='terminal exhausted')].index
        sub_t = df_nozeros.loc[cell_index].mean()* d[2]

        cell_index = metadata[(metadata['condition']==gene) & (metadata['state']=='cycling')].index
        sub_c = df_nozeros.loc[cell_index].mean()* d[3]

        f = (sub_p + sub_e + sub_t + sub_c) / (np.sum(d[0:-1]))
        temp.loc[gene] = f

    validation_test_df.loc[key] = temp.mean()

validation_df = validation_test_df.loc[['Aqr', 'Bach2', 'Bhlhe40']]
#validation_df = validation_test_df.loc['Bach2']

test_df = validation_test_df.loc[['Ets1', 'Fosb', 'Mafk', 'Stat3']]
#test_df = validation_test_df.loc['Mafk']


# %%
validation_df

# %%
test_df

# %%
test_dl = learn.dls.test_dl(validation_df, num_workers=0, shuffle=False)
preds = learn.get_preds(dl=test_dl)

clas_1 = preds[0][0]
clas_2 = preds[0][1]
clas_3 = preds[0][2]

validation_list = [clas_1.tolist(), clas_2.tolist(), clas_3.tolist()]
validation_set = pd.DataFrame(validation_list, columns=['a_i','b_i','c_i','d_i','e_i'])
validation_set.insert(0, 'genes', validation_df.index)
validation_set = validation_set.set_index('genes')
sums = validation_set.sum(axis=1)
for index in sums.index:
    sum = sums[index]
    validation_set.loc[index] = validation_set.loc[index]/sum
validation_set.to_csv('/home/kokyriakidis/Downloads/validation_output.csv')

# %%
test_dl = learn.dls.test_dl(test_df, num_workers=0, shuffle=False)
preds = learn.get_preds(dl=test_dl)

clas_1 = preds[0][0]
clas_2 = preds[0][1]
clas_3 = preds[0][2]
clas_4 = preds[0][3]

validation_list = [clas_1.tolist(), clas_2.tolist(), clas_3.tolist(), clas_4.tolist()]
test_set = pd.DataFrame(validation_list, columns=['a_i','b_i','c_i','d_i','e_i'])
test_set.insert(0, 'genes', test_df.index)
test_set = test_set.set_index('genes')
sums = test_set.sum(axis=1)
for index in sums.index:
    sum = sums[index]
    test_set.loc[index] = test_set.loc[index]/sum
test_set.to_csv('/home/kokyriakidis/Downloads/test_output.csv')


