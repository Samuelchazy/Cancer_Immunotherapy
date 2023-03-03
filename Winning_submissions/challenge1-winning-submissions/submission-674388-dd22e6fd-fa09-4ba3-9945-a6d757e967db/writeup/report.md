Predictions are made with random forest. 14 features are extracted from the sparse matrix for each of the 64 training genes and 7 testing genes. Features are
1. ratio of nonzero expression levels in unperturbed cells
2. mean of expression levels in unperturbed cells
3. mean of nonzero expression levels in unperturbed cells
4. standard deviation of expression levels in unperturbed cells
5. standard deviation of nonzero expression levels in unperturbed cells
6. skewness of expression levels in unperturbed cells
7. skewness of nonzero expression levels in unperturbed cells
8.-14. the same as above, but for perturbed cells

The output is considered as an element of R^5, so a single random forest is trained (NOT a separate RF for every state). 