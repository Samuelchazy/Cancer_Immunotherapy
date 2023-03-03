import sys
import numpy as np
import scanpy as sc

# usage: python3 calc_score.py <training file path>

states = ["progenitor", "effector", "terminal exhausted", "cycling", "other"]

# these parameters should be changed to desired values
Q = [0.5, 0.3, 0.0, 0.2, 0.0] # desired cell state proportion vector
w = [1/0.0675, 1/0.2097, 1/0.3134, 1/0.3921, 0] # use same weights as challenge 2 part B

#########################################################################################
# input: 
#   s0 - proportion vector for the unperturbed samples
#   si - proportion vector the gene i perturbed samples
#   ni - the number of gene i perturbed samples
# output: the score

def calc_score(s0, si, ni):
  if ni <= 0: print("No samples!"); return -1
  num = sum(w[j]*abs(Q[j] - si[j]) for j in range(5))
  den = sum(w[j]*abs(Q[j] - s0[j]) for j in range(5))
  if abs(den) < 0.0001: print("Unperturbed statistic matches desired state!"); return -1
  return (1 - num/den) * (1 - 1/ni**0.5)

#########################################################################################

# read the training file
print('Reading "sc_training.h5ad"')
adata = sc.read_h5ad(sys.argv[1]) 
conditions = adata.obs["condition"].unique().to_numpy()
counts = adata.obs['condition'].value_counts()

# n = the number of gene i perturbed samples
n = [counts[condition] for condition in conditions]

# s = the state proportion vectors for each gene in conditions
s = []
df = adata.obs[['condition', 'state']]
props = df.groupby('condition', as_index=False).value_counts(normalize=True)
for condition in conditions:
  condProps = props[props['condition'] == condition]  
  s.append([condProps[condProps['state'] == state].values[0][2] for state in states])
s0 = s[np.where(conditions == "Unperturbed")[0][0]]

# output scores for each gene in conditions
print("       gene:    score    n                                             s")
for i in range(len(conditions)):
  score = calc_score(s0, s[i], n[i])
  print("%11s:%9.5f%5d [%0.5f, %0.5f, %0.5f, %0.5f, %0.5f]" \
        %(conditions[i], score, n[i], s[i][0], s[i][1], s[i][2], s[i][3], s[i][4]))
