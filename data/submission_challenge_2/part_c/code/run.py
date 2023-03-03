#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!pip3 install scanpy


# # Challenge 1

# In[2]:


import re
import pandas as pd
import scanpy as sc
import os
import numpy as np
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings('ignore')


# # Read observation data

# In[3]:


os.getcwd()
os.chdir('data')


# ### Get observation data

# In[4]:


print('reading observation data...')
df_obs = sc.read_h5ad('./sc_training.h5ad')
df_obs = df_obs.obs
df_obs


# ### Unperturbed data

# In[5]:


# unperturbed_data = df_obs[df_obs['condition'] == 'Unperturbed']
# unperturbed_data = unperturbed_data[['gRNA_maxID','state']]
# unperturbed_data = unperturbed_data.groupby(by=['gRNA_maxID'],as_index=False).value_counts(normalize=True)
#
# unperturbed_data


# ### Get cells data

# In[6]:


print('reading cells adata.X.toarray()...')
adata_df = sc.read_h5ad('./sc_training.h5ad')
adata = pd.DataFrame(adata_df.X.toarray())
adata.columns = adata_df.var_names
adata.index = adata_df.obs_names
adata


# ### Get clone information file

# In[7]:


clone_info = pd.read_csv('clone_information.csv')
clone_info.index = clone_info['Unnamed: 0']
clone_info = clone_info.drop(['Unnamed: 0','gRNA_maxID'],axis=1)
clone_info = clone_info.rename_axis('index')
clone_info


# ### Merge all the data

# In[8]:


print('merging observation data with cells data on the same index...')
df_partial_data = df_obs.merge(adata,left_index=True,right_index=True)
df_full_data = df_partial_data.merge(clone_info,left_index=True,right_index=True)

df_full_data


# ### Perturbed data

# In[9]:


print('dropping unperturbed data...')
perturbed_data = df_full_data[df_full_data['condition'] != 'Unperturbed']
perturbed_data['condition'] = perturbed_data['condition'].cat.remove_categories('Unperturbed')

perturbed_data


# In[10]:


# print('grouping the data by gRNA_maxID and getting the proportions...')
# perturbed_data_prop = perturbed_data.iloc[:,0:4].groupby(by=['gRNA_maxID'],as_index=False).value_counts(normalize=True)
# perturbed_data_prop


# In[11]:


# print('dropping perturbations with < 2...')
# perturbed_data_prop_filtered = perturbed_data_prop[perturbed_data_prop['proportion'] > 2]
# perturbed_data_prop_filtered = perturbed_data_prop_filtered.drop('proportion',axis=1).reset_index(drop=True)
# perturbed_data_prop_filtered


# In[12]:


ax,fig = plt.subplots(figsize=(20,3),dpi=300)
ax = sns.countplot(x=perturbed_data['state'],color='green',alpha=0.3)
plt.grid()


# # Check null values

# In[13]:


def check_nans(df):
    nans_data = df.isna().sum().reset_index()
    nans_data = nans_data[nans_data[0] > 0]
    return nans_data


# In[14]:


check_nans(perturbed_data)


# In[15]:


perturbed_data['umi'] = perturbed_data['umi'].fillna('mode')
perturbed_data_final = perturbed_data


# In[16]:


perturbed_data_final['umi']


# # Split the data into Train & validate

# In[17]:


print('splitting the data into train, validate, & test sets...')
train_data = perturbed_data_final.drop(['Aqr','Bach2','Bhlhe40','Ets1','Fosb','Mafk','Stat3'],axis=1)
validate_data = perturbed_data_final.iloc[:,0:4]
validate_data['gRNA_bam'] = perturbed_data_final['gRNA_bam']
validate_data['umi'] = perturbed_data_final['umi']
validate_data = pd.concat([validate_data,perturbed_data_final[['Aqr','Bach2','Bhlhe40']]],axis=1)


# In[18]:


validate_data


# ### Encoding categorical columns

# In[19]:


from sklearn.preprocessing import LabelEncoder

def encode_data(df_train,df_validate):

    encoder_A = LabelEncoder()
    encoder_B = LabelEncoder()
    encoder_C = LabelEncoder()
    encoder_D = LabelEncoder()
    encoder_E = LabelEncoder()
    encoder_F = LabelEncoder()

    print('encoding train categorical variables...')
    df_train['gRNA_maxID'] = encoder_A.fit_transform(df_train['gRNA_maxID'])
    df_train['state'] = encoder_B.fit_transform(df_train['state'])
    df_train['condition'] = encoder_C.fit_transform(df_train['condition'])
    df_train['lane'] = encoder_D.fit_transform(df_train['lane'])
    df_train['gRNA_bam'] = encoder_E.fit_transform(df_train['gRNA_bam'])
    df_train['umi'] = encoder_F.fit_transform(df_train['umi'])

    print('encoding validate categorical variables...')
    df_validate['gRNA_maxID'] = encoder_A.transform(df_validate['gRNA_maxID'])
    df_validate['state'] = encoder_B.transform(df_validate['state'])
    df_validate['condition'] = encoder_C.transform(df_validate['condition'])
    df_validate['lane'] = encoder_D.transform(df_validate['lane'])
    df_validate['gRNA_bam'] = encoder_E.fit_transform(df_validate['gRNA_bam'])
    df_validate['umi'] = encoder_F.fit_transform(df_validate['umi'])

    return df_train,df_validate

train_data_labeled,validate_data_labeled = encode_data(train_data,validate_data)
train_data_labeled


# In[20]:


validate_data_labeled


# ### Normalize data

# In[21]:


def normalize_data(df_train,df_validate):
    print('normalizing training data...')
    df_train_normalized = df_train.copy()
    df_train_normalized.iloc[:,4:-2] = (df_train_normalized.iloc[:,4:-2] / df_train_normalized.iloc[:,4:-2].sum(axis=0)) * 10000
    df_train_normalized.iloc[:,4:-2] = np.log1p(df_train_normalized.iloc[:,4:-2])

    print('normalizing validate data...')
    df_validate_normalized = df_validate.copy()
    df_validate_normalized.iloc[:,6:] = (df_validate_normalized.iloc[:,6:] / df_validate_normalized.iloc[:,6:].sum(axis=0)) * 10000
    df_validate_normalized.iloc[:,6:] = np.log1p(df_validate_normalized.iloc[:,6:])

    return df_train_normalized,df_validate_normalized

df_train_data_normalized,df_validate_data_normalized = normalize_data(train_data_labeled,validate_data_labeled)
df_train_data_normalized


# In[22]:


nans_train = check_nans(df_train_data_normalized)
nans_train


# In[23]:


df_train_data_normalized = df_train_data_normalized.fillna(0)
df_train_data_normalized


# In[24]:


nans_train = check_nans(df_train_data_normalized)
nans_train


# ### Dimensionality reduction: applying PCA

# In[25]:


from sklearn.decomposition import KernelPCA

print('applying PCA of 1 component to the data...')
pca = KernelPCA(n_components=1,kernel='rbf')

left_side = df_train_data_normalized.iloc[:,:4]
left_side['gRNA_bam'] = df_train_data_normalized['gRNA_bam']
left_side['umi'] = df_train_data_normalized['umi']
left_side = left_side.reset_index(drop=True)
right_side = df_train_data_normalized.iloc[:,4:-2]
right_side = right_side.reset_index(drop=True)

right_side_pca = pd.DataFrame(pca.fit_transform(right_side))
df_train_data_normalized_PCA = pd.concat([left_side,right_side_pca],axis=1)
df_train_data_normalized_PCA


# # Split the data into X_train,X_valid,y

# In[26]:


from sklearn.model_selection import train_test_split

print('splitting the data into X & y...')
X_train = df_train_data_normalized_PCA.drop('state',axis=1).reset_index(drop=True)
X_train.columns = X_train.columns.astype(str)
y = df_train_data_normalized_PCA['state'].reset_index(drop=True)
X_valid = df_validate_data_normalized.drop('state',axis=1)

X_train


# In[27]:


X_valid


# # Balance the target variable

# In[28]:


from imblearn.over_sampling import SMOTE

smote = SMOTE(random_state=42)
X_train_res,y_res = smote.fit_resample(X_train,y)


# In[29]:


ax,fig = plt.subplots(figsize=(20,3),dpi=300)
ax = sns.countplot(x=y_res,color='green',alpha=0.3)
plt.grid()


# ### Use KMeans to determine the optimal number of clusters

# In[30]:


from sklearn.cluster import KMeans

k_values = range(1,20)
wcss = []

for k in k_values:
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(X_train_res)
    score = kmeans.inertia_
    wcss.append(score)

plt.plot(k_values, wcss)
plt.xlabel('Number of clusters')
plt.ylabel('WCSS')
plt.show()


# In[31]:


# import umap
# import seaborn as sns
# import random
#
# k = 4
#
# reducer = umap.UMAP(min_dist=0.3,n_neighbors=k,verbose=2)
# reduced_data_train = reducer.fit_transform(X_train_res,y_res)
#
# labels = y_res
# labels = labels.replace({0:'cycling',1:'effector',3:'other',4:'progenitor',2:'terminal exhausted'})
#
# ax,fig = plt.subplots(figsize=(3.5,3.5),dpi=150)
# sns.scatterplot(x=reduced_data_train[:,0],y=reduced_data_train[:,1],hue=labels,palette='magma')
# plt.legend(loc='center left',bbox_to_anchor=(1,0.5),fontsize=6)
# plt.xticks(fontsize=6)
# plt.yticks(fontsize=6)
# plt.grid();


# ### Build a K-nearest neigbors model using euclidean distance to cluster the data

# In[32]:


from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import f1_score

validate_vars = X_valid.iloc[:,5:].columns

def fit_knn_model(X_train,y,X_valid,vars):

    final_data_frame = pd.DataFrame(np.zeros([1,5]),columns=['a_i','b_i','c_i','d_i','e_i'])
    for var in vars:
        print(f'fitting, predicting, & adding the clusters of a KNeighborsClassifier with 4 neighbors for {var}...')
        knn = KNeighborsClassifier(n_neighbors=k,metric='euclidean')
        knn_model = knn.fit(X_train,y)

        X_valid_df = X_valid[['gRNA_maxID','condition','lane','gRNA_bam','umi',var]]
        X_valid_df.columns = ['gRNA_maxID','condition','lane','gRNA_bam','umi','0']
        y_pred_knn_validate = knn_model.predict(X_valid_df)
        X_valid_df['cluster'] = y_pred_knn_validate

        X_valid_proportion = X_valid_df.groupby('cluster',as_index=False)['0'].sum()
        X_valid_proportion = X_valid_proportion['0']/X_valid_proportion['0'].sum()
        X_valid_proportion = pd.DataFrame(X_valid_proportion).T
        X_valid_proportion.columns = ['a_i','b_i','c_i','d_i','e_i']
        final_data_frame = pd.concat([final_data_frame,X_valid_proportion],axis=0)

    final_data_frame = final_data_frame.drop(0,axis=0)
    final_data_frame.index = vars
    final_data_frame = pd.DataFrame(round(final_data_frame,9))
    final_data_frame = final_data_frame.rename_axis('gene')
    return final_data_frame

validate_data_predicted = fit_knn_model(X_train_res,y_res,X_valid,validate_vars)
validate_data_predicted


# In[33]:


validate_data_predicted.to_csv('submission_challenge_1/solution/validation_output.csv',index=True)


# ### Get the ground truth

# In[34]:


def get_ground_truth(df,vars):
    final_data_frame = pd.DataFrame(np.zeros([1,5]),columns=['a_i','b_i','c_i','d_i','e_i'])
    for var in vars:
        ground_truth = df.groupby('state',as_index=False)[var].sum()
        ground_truth = ground_truth[var]/ground_truth[var].sum()
        ground_truth = pd.DataFrame(ground_truth).T
        ground_truth.columns = ['a_i','b_i','c_i','d_i','e_i']
        final_data_frame = pd.concat([final_data_frame,ground_truth],axis=0)

    final_data_frame = final_data_frame.drop(0,axis=0)
    final_data_frame.index = vars
    final_data_frame = round(final_data_frame,9)
    return final_data_frame

ground_truth = get_ground_truth(validate_data,validate_vars)
ground_truth


# ### Compare the prediction to the ground truth

# In[35]:


diff = validate_data_predicted - ground_truth
diff = diff.apply(lambda x: round(x * 100,9))
diff = diff.apply(lambda x: x.astype(str)+' %' )
diff = diff.rename_axis(['diff_in_%'])
diff


# # Predict test data

# In[46]:


print('splitting the data into train,validate, & test sets...')
X = perturbed_data_final.drop(['Ets1','Fosb','Mafk','Stat3'],axis=1)
test_data = pd.concat([perturbed_data_final.iloc[:,0:4],perturbed_data_final[['Ets1','Fosb','Mafk','Stat3']]],axis=1)
test_data['gRNA_bam'] = perturbed_data_final['gRNA_bam']
test_data['umi'] = perturbed_data_final['umi']
X


# In[47]:


check_nans(test_data)


# In[48]:


test_data


# ### Encoding categorical columns

# In[49]:


from sklearn.preprocessing import LabelEncoder

def encode_data(X_df,df_test):

    encoder_A = LabelEncoder()
    encoder_B = LabelEncoder()
    encoder_C = LabelEncoder()
    encoder_D = LabelEncoder()
    encoder_E = LabelEncoder()
    encoder_F = LabelEncoder()

    print('encoding train categorical variables...')
    X_df['gRNA_maxID'] = encoder_A.fit_transform(X_df['gRNA_maxID'])
    X_df['state'] = encoder_B.fit_transform(X_df['state'])
    X_df['condition'] = encoder_C.fit_transform(X_df['condition'])
    X_df['lane'] = encoder_D.fit_transform(X_df['lane'])
    X_df['gRNA_bam'] = encoder_E.fit_transform(X_df['gRNA_bam'])
    X_df['umi'] = encoder_F.fit_transform(X_df['umi'])

    print('encoding test categorical variables...')
    df_test['gRNA_maxID'] = encoder_A.transform(df_test['gRNA_maxID'])
    df_test['state'] = encoder_B.transform(df_test['state'])
    df_test['condition'] = encoder_C.transform(df_test['condition'])
    df_test['lane'] = encoder_D.transform(df_test['lane'])
    df_test['gRNA_bam'] = encoder_E.fit_transform(df_test['gRNA_bam'])
    df_test['umi'] = encoder_F.fit_transform(df_test['umi'])

    return X_df,df_test

X_data_labeled,test_data_labeled = encode_data(X,test_data)
X_data_labeled


# In[50]:


test_data_labeled


# ### Normalize training data

# In[51]:


def normalize_data(df_X,df_test):
    print('normalizing training data...')
    df_train_normalized = df_X.copy()
    df_train_normalized.iloc[:,4:-2] = (df_train_normalized.iloc[:,4:-2] / df_train_normalized.iloc[:,4:-2].sum(axis=0)) * 10000
    df_train_normalized.iloc[:,4:-2] = np.log1p(df_train_normalized.iloc[:,4:-2])

    print('normalizing test data...')
    df_test_normalized = df_test.copy()
    df_test_normalized.iloc[:,4:-2] = (df_test_normalized.iloc[:,4:-2] / df_test_normalized.iloc[:,4:-2].sum(axis=0)) * 10000
    df_test_normalized.iloc[:,4:-2] = np.log1p(df_test_normalized.iloc[:,4:-2])

    return  df_train_normalized,df_test_normalized

df_X_data_normalized,df_test_normalized = normalize_data(X_data_labeled,test_data_labeled)
df_X_data_normalized


# In[52]:


df_test_normalized


# In[53]:


nans_X = check_nans(df_X_data_normalized)
nans_X


# In[54]:


df_X_data_normalized = df_X_data_normalized.fillna(0)
df_X_data_normalized


# In[55]:


check_nans(df_test_normalized)


# In[56]:


df_X_data_normalized = df_X_data_normalized.drop(nans_X['index'],axis=1)
df_X_data_normalized


# ### Dimensionality reduction: applying PCA

# In[57]:


from sklearn.decomposition import KernelPCA

print('applying PCA of 1 component to the data...')
pca = KernelPCA(n_components=1,kernel='rbf')

left_side = df_test_normalized.iloc[:,:4]
left_side['gRNA_bam'] = df_test_normalized['gRNA_bam']
left_side['umi'] = df_test_normalized['umi']
left_side = left_side.reset_index(drop=True)
right_side = df_test_normalized.iloc[:,4:-2]
right_side = right_side.reset_index(drop=True)
right_side_pca = pd.DataFrame(pca.fit_transform(right_side))
df_X_data_normalized_PCA = pd.concat([left_side,right_side_pca],axis=1)
df_X_data_normalized_PCA


# ### Split the data into X,X_test,y

# In[58]:


print('splitting the data into X & y...')
X_train = df_X_data_normalized_PCA.drop('state',axis=1).reset_index(drop=True)
X_train.columns = X_train.columns.astype(str)
y_train = df_X_data_normalized_PCA['state'].reset_index(drop=True)
X_test = df_test_normalized.drop('state',axis=1)

X_train


# In[59]:


from imblearn.over_sampling import SMOTE

smote = SMOTE(random_state=42)
X_train_res,y_train_res = smote.fit_resample(X_train,y_train)


# In[60]:


ax, fig = plt.subplots(figsize=(20, 3), dpi=300)
ax = sns.countplot(x=y_train_res, color='green', alpha=0.3)
plt.grid()


# In[61]:


X_train_res


# In[62]:


test_data


# In[63]:


test_vars = test_data.iloc[:,4:-2].columns

def fit_knn_model(X,y,X_test,vars):

    final_data_frame = pd.DataFrame(np.zeros([1,5]),columns=['a_i','b_i','c_i','d_i','e_i'])
    for var in vars:
        print(f'fitting, predicting, & adding the clusters of a KNeighborsClassifier with 4 neighbors for {var}...')
        knn = KNeighborsClassifier(n_neighbors=k,metric='euclidean')
        knn_model = knn.fit(X,y)

        X_test_df = X_test[['gRNA_maxID','condition','lane','gRNA_bam','umi',var]]
        X_test_df.columns = ['gRNA_maxID','condition','lane','gRNA_bam','umi','0']
        y_pred_knn_test = knn_model.predict(X_test_df)
        X_test_df['cluster'] = y_pred_knn_test

        X_test_proportion = X_test_df.groupby('cluster',as_index=False)['0'].sum()
        X_test_proportion = X_test_proportion['0']/X_test_proportion['0'].sum()
        X_test_proportion = pd.DataFrame(X_test_proportion).T
        X_test_proportion.columns = ['a_i','b_i','c_i','d_i','e_i']
        final_data_frame = pd.concat([final_data_frame,X_test_proportion],axis=0)

    final_data_frame = final_data_frame.drop(0,axis=0)
    final_data_frame.index = vars
    final_data_frame = pd.DataFrame(round(final_data_frame,9))
    final_data_frame = final_data_frame.rename_axis('gene')
    return final_data_frame

test_data_predicted = fit_knn_model(X_train_res,y_train_res,test_data,test_vars)
test_data_predicted


# In[64]:


test_data_predicted.to_csv('submission_challenge_1/solution/test_output.csv',index=True)


# ### get the ground truth

# In[65]:


ground_truth_test = get_ground_truth(test_data,test_vars)
ground_truth_test


# ### Compare the prediction to the ground truth

# In[66]:


diff_test = test_data_predicted - ground_truth_test
diff_test = diff_test.apply(lambda x: round(x * 100,9))
diff_test = diff_test.apply(lambda x: x.astype(str)+' %' )
diff_test = diff_test.rename_axis(['diff_in_%'])
diff_test


# ### Read guide abundance f ile

# In[67]:


# guide_abundance = pd.read_csv('guide_abundance.csv')
# guide_abundance.columns = ['guide_ID','plasmid_pool','perturbseq']
# guide_abundance


# In[68]:


df_full_data


# ### Read metadata

# In[69]:


# scRNA_ATAC = sc.read_10x_h5('./scRNA_ATAC.h5')


# In[70]:


# scRNA = scRNA_ATAC.copy()
# scRNA = scRNA.var.reset_index(drop=False)
# scRNA = scRNA[['index','gene_ids']]
# scRNA.columns = ['condition','gene_ids']
# scRNA = scRNA.sort_values(by=('condition')).reset_index(drop=True)
# scRNA


# # Challenge 2

# In[71]:


df_total_data = perturbed_data_final.copy()
df_total_data


# ### Encoding categorical columns

# In[72]:


def encode_data(df):

    encoder_A = LabelEncoder()
    encoder_B = LabelEncoder()
    encoder_C = LabelEncoder()
    encoder_D = LabelEncoder()
    encoder_E = LabelEncoder()
    encoder_F = LabelEncoder()

    print('encoding categorical variables...')
    df['gRNA_maxID'] = encoder_A.fit_transform(df['gRNA_maxID'])
    df['state'] = encoder_B.fit_transform(df['state'])
    df['condition'] = encoder_C.fit_transform(df['condition'])
    df['lane'] = encoder_D.fit_transform(df['lane'])
    df['gRNA_bam'] = encoder_E.fit_transform(df['gRNA_bam'])
    df['umi'] = encoder_F.fit_transform(df['umi'])

    return df

df_total_data_labeled = encode_data(df_total_data)
df_total_data_labeled


# ### Normalizing data

# In[73]:


def normalize_data(df):
    print('normalizing data...')
    df_train_normalized = df.copy()
    df_train_normalized.iloc[:,4:-2] = (df_train_normalized.iloc[:,4:-2] / df_train_normalized.iloc[:,4:-2].sum(axis=0)) * 10000
    df_train_normalized.iloc[:,4:-2] = np.log1p(df_train_normalized.iloc[:,4:-2])

    return df_train_normalized

df_total_data_labeled_normalized = normalize_data(df_total_data_labeled)
df_total_data_labeled_normalized


# ### Check null values

# In[77]:


nans_total_data = check_nans(df_total_data_labeled_normalized)
nans_total_data


# In[78]:


df_total_data_labeled_normalized = df_total_data_labeled_normalized.fillna(0)
df_total_data_labeled_normalized


# ### Dimensionality reduction: applying PCA

# In[79]:


from sklearn.decomposition import KernelPCA

print('applying PCA of 1 component to the data...')
pca = KernelPCA(n_components=1,kernel='rbf')

left_side = df_total_data_labeled_normalized.iloc[:,:4]
left_side['gRNA_bam'] = df_total_data_labeled_normalized['gRNA_bam']
left_side['umi'] = df_total_data_labeled_normalized['umi']
left_side = left_side.reset_index(drop=True)
right_side = df_total_data_labeled_normalized.iloc[:,4:-2]
right_side = right_side.reset_index(drop=True)
right_side_pca = pd.DataFrame(pca.fit_transform(right_side))
df_total_data_labeled_normalized_PCA = pd.concat([left_side,right_side_pca],axis=1)
df_total_data_labeled_normalized_PCA


# # Split the Data into X & y

# In[80]:


print('splitting the data into X & y...')
X_final = df_total_data_labeled_normalized_PCA.drop('state',axis=1).reset_index(drop=True)
X_final.columns = X_final.columns.astype(str)
y_final = df_total_data_labeled_normalized_PCA['state'].reset_index(drop=True)

X_final


# In[81]:


from imblearn.over_sampling import SMOTE

smote = SMOTE(random_state=42)
X_final_res,y_final_res = smote.fit_resample(X_final,y_final)

X_final_res


# In[82]:


check_nans(X_final_res)


# In[83]:


ax, fig = plt.subplots(figsize=(20, 3), dpi=300)
ax = sns.countplot(x=y_final_res, color='green', alpha=0.3)
plt.grid()


# ### Filter the genes

# In[84]:


df_genes = df_total_data_labeled_normalized.iloc[:,4:-2]
df_genes


# # Predict the final state of all genes (around 1.5 hours)

# In[85]:


def fit_knn_model(X,y,df_genes):

    final_data_frame = pd.DataFrame(np.zeros([1,5]),columns=['a_i','b_i','c_i','d_i','e_i'])
    X.columns = X.columns.astype(str)
    for gene in tqdm(df_genes.columns):
        knn = KNeighborsClassifier(n_neighbors=k,metric='euclidean')
        knn_model = knn.fit(X,y)

        X_df = df_total_data_labeled_normalized[['gRNA_maxID','condition','lane','gRNA_bam','umi',gene]]
        X_df.columns = ['gRNA_maxID','condition','lane','gRNA_bam','umi','0']
        y_pred_knn_test = knn_model.predict(X_df)
        X_df['cluster'] = y_pred_knn_test

        X_proportion = X_df.groupby('cluster',as_index=False)['0'].sum()
        X_proportion = X_proportion['0']/X_proportion['0'].sum()
        X_proportion = pd.DataFrame(X_proportion).T
        X_proportion.columns = ['a_i','b_i','c_i','d_i','e_i']
        final_data_frame = pd.concat([final_data_frame,X_proportion],axis=0)

    final_data_frame = final_data_frame.drop(0,axis=0)
    final_data_frame.index = df_genes.columns
    final_data_frame = pd.DataFrame(round(final_data_frame,9))
    final_data_frame = final_data_frame.rename_axis('gene')
    return final_data_frame

final_data_predicted = fit_knn_model(X_final_res,y_final_res,df_genes)
final_data_predicted.to_csv('final_data_predicted.csv',index=True)
final_data_predicted


# In[86]:


final_data = pd.read_csv('final_data_predicted.csv')
final_data


# In[87]:


check_nans(final_data)


# In[88]:


final_data = final_data.fillna(0)
final_data


# # Challenge part_a

# In[89]:


final_data_a_ordered = final_data.copy()

final_data_a_ordered = final_data_a_ordered.sort_values(by='a_i',ascending=False)
final_data_a_ordered['cycling_constraint'] = np.where(final_data_a_ordered['d_i'].values > 0.05,1,0)
final_data_a_ordered = final_data_a_ordered.rename(columns={'a_i':'objective'})
final_data_a_ordered.index = final_data_a_ordered['gene']
final_data_a_ordered = final_data_a_ordered[['objective','cycling_constraint']]

final_data_a_ordered


# ### Submission part_a

# In[90]:


final_data_a_ordered.to_csv('submission_challenge_2/part_a/part_a_output.csv',index=True)


# # Challenge part_b

# In[91]:


final_data_b_ordered = final_data.copy()

final_data_b_ordered['objective'] = (final_data_b_ordered['a_i'] / 0.0675) + (final_data_b_ordered['b_i'] / 0.2097) - (final_data_b_ordered['c_i'] / 0.3134) + (final_data_b_ordered['d_i'] / 0.3921)
final_data_b_ordered['objective'] = round(final_data_b_ordered['objective'],9)
final_data_b_ordered = final_data_b_ordered.sort_values(by='objective',ascending=False)
final_data_b_ordered['cycling_constraint'] = np.where(final_data_b_ordered['d_i'].values > 0.05,1,0)
final_data_b_ordered.index = final_data_b_ordered['gene']
final_data_b_ordered = final_data_b_ordered[['objective','cycling_constraint']]

final_data_b_ordered


# ### Submission part_b

# In[92]:


final_data_b_ordered.to_csv('submission_challenge_2/part_b/part_b_output.csv',index=True)


# # Challenge part_c

# ### Submission part_c

# In[93]:


final_data.to_csv('submission_challenge_2/part_c/part_c_output.csv',index=True)


# In[ ]:




