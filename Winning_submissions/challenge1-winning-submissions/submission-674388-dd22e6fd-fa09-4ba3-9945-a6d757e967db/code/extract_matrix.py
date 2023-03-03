import scanpy as sc
import scipy
from warnings import filterwarnings
filterwarnings('ignore')

sc.settings.verbosity = 3
adata = sc.read_h5ad('./sc_training.h5ad')

xf = scipy.sparse.find(adata.X)

f = open('xf.csv', 'w')

N = 80497219

for i in range(N):
    if (i % 10000 == 0): print(i, end='\r')
    print(xf[0][i], xf[1][i], xf[2][i], file=f, sep=',')

f.close()