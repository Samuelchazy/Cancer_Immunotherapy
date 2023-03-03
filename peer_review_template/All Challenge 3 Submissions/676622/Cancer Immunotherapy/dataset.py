
#  By Miftah Rangkuti

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from csv import writer

with h5py.File('sc_training.h5','r') as hdf:
    ls = list (hdf.keys())
    # print('List of datasets in this file: \n',ls)
    #list of datasets in this file:
    #['dataset1','dataset2']

    data = hdf.get('dataset1')

# listdata = list(data)
listTox2, listTpt1, listTcf7, listDvf2, listZeb2, listSox4 = []
nontarget = []

# obs1Tox2 = list.find('CAACGGT-1')
# obs2Tox2 = list.find('TCGTTTG-1')
# obs3Tox2 = list.find('AGGTTCA-1')

# obs1Tpt1 = list.find('TACGACG-1')

# obs1Tcf7 = list.find('TGAAGAG-1')

# obs1Dvf2 = list.find('AGTGTTG-1')
# obs1Dvf2 = list.find('TTTCCTC-1')
    
# obs1Zeb2 = list.find('CTCGCAT-1')

# obs1Sox4 = list.find('TTCAACT-1')
# end so on....it should be make a Datafile like JSON

while True:   
    listdata = list(data)
    if (listdata == list.find('CAACGGT-1')) or (listdata ==  list.find('CAACGGT-1')) or (listdata==list.find('AGGTTCA-1')) :
        listTox2.append(listdata)
    else:
        if (listdata == list.find('TACGACG-1')):
            listTpt1.append(listdata)
        else:
            if (listdata == list.find('TGAAGAG-1')) :
                listTcf7.append(listdata)
            else: 
                if (listdata == list.find('AGTGTTG-1')) or (listdata == list.find('TTTCCTC-1')):
                    listDvf2.append(listdata)
                else:
                    if (listdata == list.find('CTCGCAT-1')):
                        listZeb2.append(listdata)
                    else:
                        if (listdata == list.find('TTCAACT-1')):
                            listSox4.append(listdata) 
                        else:
                            nontarget.append(listdata)
       
        break


    f=h5py.File('sc_training.h5','r')
    ls = list(f.keys())
    f.close()

with open('validation_output.csv','w', encoding='utf8', newline='') as filecsv:
    thewriter = writer(filecsv)
    # header = ['aData','sgRNA_max']
    # thewriter.writerow(header)

    # savetocsv_class(listTox2)
    # savetocsv_class(listTpt1)
    # savetocsv_class(listTcf7) 
    # savetocsv_class(listDvf2)
    # savetocsv_class(listZeb2)
    # savetocsv_class(listSox4)

    for list in listTox2:
        datamatrix = np.arange(listTox2).reshape(3, 5)
        dataTox2 = list.find(datamatrix).text.replace('\n','')
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
    # make class.

# def savetocsv_class(listtype):
#     for list in listtype:
#         datamatrix = np.arange(listtype).reshape(3, 5)
#         listtype = list.find(datamatrix).text.replace('\n', '')
        # print(listtype)
        # thewriter.writerow(listtype)

