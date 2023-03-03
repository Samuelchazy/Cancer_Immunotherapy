## Report

My baseline model was a very simple model which would always predict the mean of the training targets. I wasn't able to beat this baseline while using the expression data supplied by the challenge host.

Finally, I used information about exons which can be found on the [National Library of Medicine website.](https://www.ncbi.nlm.nih.gov/gene)

The features were limited to a basic statistical description of the lengths of the exons (minimum, maximum, mean, median, first quartile and third quartile of the lengths) of each gene.

Only conditions which have at least 25 cells were selected as targets (only 55 samples).

The models were evaluated using a Leave One Out cross-validation procedure. Finally, an extra-trees regressor (from the scikit-learn library) was used as the model.