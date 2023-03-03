# writeup.md
# Cancer Immunotherapy Data Science Grand Challenge - Challenge 1

~~~
Author: Brody Langille, Jordan Trajkovski
Date: 2023-01-27
~~~

We broke down and approached the Challenge 1 problem as follows:

## Given

The following information was given about the challenge 1 problem.

* The experiment detailed in challenge 1 generated a single dataset that can be reconstructed from the provided data file
    * The .h5ad file can be read using the Python `scanpy` package, into an AnnData object
        * The AnnData.obs and anndata.X objects contain the required information and can be joined into a single dataset
    * This dataset has the shape of [28697 rows x 15079 columns]
        * 1 row per T cell/sample (28,697 total rows):
            * These T cells were injected into and recovered from a mouse with a melanoma tumor, after a 2 week period
            * Some of these T cells had specific genes "knocked out" using CRISPR prior to injection where other T cells were "unperturbed" i.e. no gene edits made prior to injection
            * Each row has 15,079 columns/values
                * 15,077 columns represent single gene expression measurements from the T cells after recovery from the mouse
                * 2 columns are labels for each T cell/sample:
                    1. Which gene was knocked out in that specific cell (using CRISPR) as a string ex. Rps6, and
                    1. What the resulting end state of that T cell was (tagged by experts) as a string ex. progenitor
    * Note: the dataset was already pre-cleaned (by others) to remove any process-related anomalies

## Problem

The problem to be solved in challenge 1 was, given only the data noted above:

* For a given list of genes to be "knocked out" - 1 at a time
* Predict the probability distribution of the T cell end states across 5 pre-defined states [progenitor, effector, terminal exhausted, cycling, other] for each of those "knocked out" genes
* Note: these genes were not already knocked out in the original experiment

## Solution

In order to predict the knockout gene effects on the T cell end state distribution we did the following steps.

1. Defined the solution as:
    * `input -> model -> output` where:
        * Input: a vector describing each "knocked out" gene in some way
        * Output: a vector describing the probability distribution of T cell end states given the (single) "knocked out" gene (vector input)
        * Model: a predictive model that allows us to predict the output, given the input, for any "knocked out" gene
1. Designed specific solutions for the input, model and output defined above:
    * `output`
        * We know the end state probability distributions for each condition/knocked out gene
        * We aggregated the single, combined input data set to a shape of [67 rows x 5 columns]
            * [67] rows:
                * 1 per "condition" i.e. "knocked out" gene + the "Unperturbed" condition
            * [5] columns:
                * 1 per "state" ex. progenitor, effector
                * Column values are the relative % of T cells/samples for that condition in that state i.e. a probability distribution across 5 total states
        * We have these `output` labels (5-dim vectors) for the 67 conditions (experimentally knocked out genes + "Unperturbed")
            * NOTE:
                * Since we do not have a corresponding input vector for "Unperturbed", since it is a not a gene, it is dropped here
                * Also, because pre-cleaning of the data (by others) resulted in no expression data available for (knocked out) genes Fzd1 and P2rx7 we also omit those 2 genes from the `output` label set (and subsequently from the matching `input` data set)
            * So, our output data set becomes shape [64 rows x 5 columns]
                * 64 genes sub-selected from the original (67) conditions
                * With 5 dimensions each (T cell end state probabilities)  
    * `input`
        * We have enough data to create vectors to describe each of the 15,077 genes in (67 dim - 1 per condition/knocked out gene) vectors
        * We aggregated the single, combined input data set to a shape of [28,697 rows x 15,077 columns]
            * [28,697] rows:
                * 1 per sample (T cell) where either: 1 gene was "knocked out" or the T cell was "Unperturbed"
            * [15,077] columns:
                * 1 per measured gene (expression)
                * Column values are the (normalized) gene expression values for each measured gene (for each row/T cell sample)
        * Since we have the `output` (labels) for the 64 _remaining_ conditions (experimentally knocked out genes) we want to use the corresponding inputs describing those knocked out genes that would generate those outputs
            * So, our input data set becomes shape [64 rows x 28,697 columns]
            * Note: although the "Unperturbed" condition is not included in the "rows" in the training data, it is still retained in the input "column" data, since we have expression data for that condition across all genes
    * `model`
        * For the model, we decided to build and train an artificial neural network (ANN) using the `input` and `output` vectors described above
            * In order to obtain the required 5 (T cell end) state probability distribution we used a Softmax activation function on the output layer of the ANN
        * The (64 x 67) `input` vectors noted above were used as the input/feature vector training data for the model
        * The (64 x 5) `output` vectors noted above were used as the output/label training data for the model
        * We iterated on various model configurations, topologies and hyper parameters to optimize our model's loss and accuracy
        * To generate our output, we:
            * Selected the (7) test and validation gene vectors from the aggregate (mean gene expression) set used for the model's input/feature vector training data
            * These vectors were then fed, one at a time, into the model to generate T cell end state probability distribution predictions
            * This full loop, executed a configurable number of times: 
                * Shuffle and segregate all input data into a 'test' set and a 'training' set, to attempt to mitigate overfitting
                * Train the model on the 'training' set ('test' samples are held out)
                * For each test and validation gene, select the input vector (28,697 expression measurements), generate a prediction with the current model, and add the result to the aggregate set for that gene
            * The final result (T cell end state probability distribution) for each of the (7) knocked out genes is an individual result that was chosen as the closest prediction to the median of the set of predictions for that gene
            * The final results were then output to their respective .csv files in the `solution/` folder as specified by the challenge
