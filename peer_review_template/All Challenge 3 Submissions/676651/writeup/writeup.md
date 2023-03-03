
# Gene Knock-Out Scoring Function for Cancer Immunotherapy

## Introduction

The goal of this write-up is to propose a statistic and a scoring function for evaluating the gene knock-out strategy in cancer immunotherapy. The desired cell state proportion vector Q represents the ideal proportions of different cell types in the treated population and is defined as (0.6, 0.25, 0.05, 0.1, 0.00), where the first component represents the proportion of progenitor cells, the second component represents the proportion of effector cells, the third component represents the proportion of terminally exhausted cells, the fourth component represents the proportion of cycling cells, and the fifth component represents the proportion of cancer cells.

## Statistic to Summarize the Gene Expression Distribution

We use the proportion of cells in each state as the statistic to summarize the gene expression distribution obtained after knocking out a gene. The proportion of cells in each state can be calculated based on the gene expression profiles of the treated population.

## Scoring Function

The score for knocking out a gene is calculated as the product of two terms: the Mahalanobis distance between the gene expression distribution of the treated population (Pi) and the gene expression distribution of the unperturbed cells (P0), and the cosine similarity between the desired cell state proportion vector Q and the proportion of cells in each state in the treated population (s(Pi)).

The Mahalanobis distance between two distributions measures the similarity between the distributions and is calculated as the square root of the sum of the squared differences between the means of the distributions, divided by the covariance matrix of the distributions. A smaller value of the Mahalanobis distance indicates a more similar distribution.

The cosine similarity between two vectors measures the cosine of the angle between the vectors and ranges from -1 to 1. A cosine similarity of 1 indicates that the vectors are identical, while a cosine similarity of -1 indicates that the vectors are orthogonal. A higher value of the cosine similarity indicates a more similar proportion of cells in each state in the treated population and the desired cell state proportion vector Q.


Here is the LaTeX formula for the score:

$$Score(P_0, Q, \hat{s}(P_i)) = (1 - D_{Mahalanobis}(P_i, P_0)) \cdot cosine\_similarity(Q, \hat{s}(P_i))$$

$$Q = (0.6, 0.25, 0.05, 0.1, 0.00)$$

Where $D_{Mahalanobis}(P_i, P_0)$ represents the Mahalanobis distance between the gene expression distribution $P_i$ obtained from knocking out gene i and the empirical gene expression distribution of the unperturbed cells $P_0$, and $cosine\_similarity(Q, \hat{s}(P_i))$ represents the cosine similarity between the desired cell state proportion vector $Q$ and the predicted statistic $\hat{s}(P_i)$. The Mahalanobis distance is calculated as follows:

$$D_{Mahalanobis}(P_i, P_0) = \sqrt{(P_i - P_0)^T \Sigma^{-1} (P_i - P_0)}$$

Where $\Sigma^{-1}$ is the inverse covariance matrix of the gene expression distribution of the unperturbed cells. $^T$ is the transpose of a matrix.

A higher value of the score indicates a better gene knock-out strategy, as it results in a more similar gene expression distribution and a more similar proportion of cells in each state in the treated population to the desired cell state proportion vector Q.

## Conclusion

In this write-up, we proposed a statistic and a scoring function for evaluating the gene knock-out strategy in cancer immunotherapy. The desired cell state proportion vector Q represents the ideal proportions of different cell types in the treated population, and the score for knocking out a gene is calculated as the product of the Mahalanobis distance between the gene expression distributions and the cosine similarity between the desired cell state proportion vector and the proportion of cells in each state in the treated population. A higher value of the score indicates a better gene knock-out strategy.