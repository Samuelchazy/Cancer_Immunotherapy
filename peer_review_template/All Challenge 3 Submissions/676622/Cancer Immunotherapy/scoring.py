import scanpy as sc
import numpy as np
import pandas as pd
import h5py
from csv import writer

adata = sc.read_h5ad('Tcell_challange.h5ad')
df = adata.obs[['gRNA_max', 'cell_state']]
props = df.groupby('gRNA_max', as_index=False).value_counts(normalize=True)

listTox2, listTpt1, listTcf7, listDvf2, listZeb2, listSox4 = []
nontarget = []

#

while True:
    listdata = list(df)
    if (listdata == list.find('CAACGGT-1')) or (listdata == list.find('CAACGGT-1')) or (listdata == list.find('AGGTTCA-1')):
        # cell1
        propsTox2 = df.value_counts(normalize=True)
        if (propsTox2 == 0.00):
            listTox2.append(listdata)

    else:
        if (listdata == list.find('TACGACG-1')):
            # cell2
            propsTpt1 = df.value_counts(normalize=True)
        if (propsTpt1 == 0.00):
            listTpt1.append(listdata)
        else:
            if (listdata == list.find('TGAAGAG-1')):
                # cell3
                propsTcf7 = df.value_counts(normalize=True)
                if (propsTcf7 == 0.00):
                    listTcf7.append(listdata)
            else:
                if (listdata == list.find('AGTGTTG-1')) or (listdata == list.find('TTTCCTC-1')):
                    # cell4
                    propsDvf2 = df.value_counts(normalize=True)
                    if (propsDvf2 == 0.00):
                        listDvf2.append(listdata)
                else:
                    if (listdata == list.find('CTCGCAT-1')):
                        # cell5
                        propsZeb2 = df.value_counts(normalize=True)
                        if (propsZeb2 == 0.00):
                            listZeb2.append(listdata)
                    else:
                        if (listdata == list.find('TTCAACT-1')):
                            # cell6
                            propsSox4 = df.value_counts(normalize=True)
                            if (propsSox4 == 0.00):
                                listSox4.append(listdata)
                        
        break

    f = h5py.File('Tcell_challange.h5ad')
    f.close()
    print(listTox2, listTpt1, listTcf7, listDvf2, listZeb2, listSox4)
