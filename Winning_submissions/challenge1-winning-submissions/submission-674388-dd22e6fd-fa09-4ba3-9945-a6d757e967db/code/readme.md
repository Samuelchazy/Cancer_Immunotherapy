DEPLOYMENT INSTRUCTIONS
=======================

- Compile *.cpp c++ source files with gcc compiler (each separately) to create *.exe
- run the consecutive 6 commands:
1. python extract_obs.py
2. python extract_var.py
3. python extract_matrix.py
4. binarize_sparse_matrix.exe
5. extract_features.exe
6. multi_rf.exe

Explanation
===========
1., 2., 3. - extracting data from sc_training.h5ad to CSV files
4. - converting sparse matrix into binary form (bindata/matrix.bin) for fast loading in C++
5. - calculation of features and ground truth (stored into training.csv and testing.csv)
6. - training random forest and predicting on test genes
