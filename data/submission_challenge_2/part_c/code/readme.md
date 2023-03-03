##### Deployment Instructions

# The original notebook was developed using Pycharm / python 3.10
# Libraries to install / import:
# scanpy, pandas, os, numpy, seaborn, tqdm, matplotlib, sklearn, imblearn, umap, random

Process:

Challenge 1:

1- Import / read all available files
2- Merge all the files into 1 file
3- Drop Unperturbed data from the condition column
4- Check the distribution of variables and especially the target variable
5- Check for null values and fill missing values
6- Split the data into Train & Validate and keep the Test dataset aside
7- Manage the Train & Validate datasets
8- Encode the categorical variables using Label Encoder
9- Normalize the data, multiply by 10,000 and apply log1p
10- Apply a kernelPCA with 1 component to all the gene columns to capture the maximum variance from these variables
11- The final dataset now has all the initial columns (gRNA_maxID, condition, lane, gRNA_bam, umi, state) + 1 combined gene variance column
12- Split the data into X_train, X_valid, and y, where y is the target feature representing the state of the gene
13- Using SMOTE, balance the target variable proportions
14- Use KMeans to determine the number of possible clusters on the X_train data
15- Build a K-nearest neighbors model on the X_train data using euclidean distance to cluster the data based on the number obtained by KMeans
16- Predict the clusters of X_validate datasets
17- Group the data by cluster, sum the values, and normalize them => the predicted vector now is the probability of the state of each gene
18- Get the ground truth which is the validate data, and compare them to the predicted data to validate the prediction accuracy
19- Predict the Test dataset by applying the same steps as before starting with step 8 by encoding the categorical variables up till step 18
20- save all the files to disk

Challenge 2:

- Repeat the previous process on all the dataset starting with step 8 by encoding the categorical variables up till step 17 by predicting all the genes probability states

Challenge part_a

- Sort the values of the predicted genes table by the progenitor state “a_i” in a descending manner
- Add a constraint column where “cycling_constraint” = 1 if the cycling state values “d_i” is bigger than 0.05, or else = 0
- save the file as part_a_output

Challenge part_b

- Add an objective function to the previous dataset with the following formula: (a_i / 0.0675) + (b_i / 0.2097) - (c_i / 0.3134) + (d_i / 0.3921)
- Sort the dataset by the objective function in a descending manner
- Keep the same cycling constraint as in part_a
- save the file as part_b_output

Challenge part_c

- This part refers to the obtained predicted dataset with all the genes and their state proportions
- save the file as part_c_output
