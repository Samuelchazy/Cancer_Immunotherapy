## CAR T-Cell Therapy

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
        if (propsTox2 >= 0.05):
            listTox2.append(listdata)

    else:
        if (listdata == list.find('TACGACG-1')):
            # cell2
            propsTpt1 = df.value_counts(normalize=True)
        if (propsTpt1 >= 0.05):
            listTpt1.append(listdata)
        else:
            if (listdata == list.find('TGAAGAG-1')):
                # cell3
                propsTcf7 = df.value_counts(normalize=True)
                if (propsTcf7 >= 0.05):
                    listTcf7.append(listdata)
            else:
                if (listdata == list.find('AGTGTTG-1')) or (listdata == list.find('TTTCCTC-1')):
                    # cell4
                    propsDvf2 = df.value_counts(normalize=True)
                    if (propsDvf2 >= 0.05):
                        listDvf2.append(listdata)
                else:
                    if (listdata == list.find('CTCGCAT-1')):
                        # cell5
                        propsZeb2 = df.value_counts(normalize=True)
                        if (propsZeb2 >= 0.05):
                            listZeb2.append(listdata)
                    else:
                        if (listdata == list.find('TTCAACT-1')):
                            # cell6
                            propsSox4 = df.value_counts(normalize=True)
                            if (propsSox4 >= 0.05):
                                listSox4.append(listdata)
                        else:
                            nontarget.append(listdata)

        break

    f = h5py.File('Tcell_challange.h5ad')
    f.close()

    with open('part_b_output.csv', 'w', encoding='utf8', newline='') as filecsv:
        thewriter = writer(filecsv)
    for list in listTox2:
        datamatrix = np.arange(listTox2).reshape(3, 5)
        dataTox2 = list.find(datamatrix).text.replace('\n', '')
        # print(dataTox2)
        thewriter.writerow(dataTox2)

    for list in listTpt1:
        datamatrix = np.arange(listTpt1).reshape(3, 5)
        listTpt1 = list.find(datamatrix).text.replace('\n', '')
        # print(listTpt1)
        thewriter.writerow(listTpt1)

    for list in listTcf7:
        datamatrix = np.arange(listTcf7).reshape(3, 5)
        listTcf7 = list.find(datamatrix).text.replace('\n', '')
        # print(listTcf7)
        thewriter.writerow(listTcf7)

    for list in listDvf2:
        datamatrix = np.arange(listDvf2).reshape(3, 5)
        listDvf2 = list.find(datamatrix).text.replace('\n', '')
        # print(listDvf2)
        thewriter.writerow(listDvf2)

    for list in listZeb2:
        datamatrix = np.arange(listZeb2).reshape(3, 5)
        listZeb2 = list.find(datamatrix).text.replace('\n', '')
        # print(listZeb2)
        thewriter.writerow(listZeb2)

    for list in listSox4:
        datamatrix = np.arange(listSox4).reshape(3, 5)
        listSox4 = list.find(datamatrix).text.replace('\n', '')
        # print(listTcf7)
        thewriter.writerow(listSox4)
