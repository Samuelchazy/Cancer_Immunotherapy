import scanpy as sc
from warnings import filterwarnings
filterwarnings('ignore')

sc.settings.verbosity = 3
adata = sc.read_h5ad('./sc_training.h5ad')
adata.var.to_csv('var.csv')