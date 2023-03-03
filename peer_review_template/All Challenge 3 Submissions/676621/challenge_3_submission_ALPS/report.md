# TEAM ALPS


## MIT CHALLENGE 3


Our approach consists of the following steps:
- We are given two cell populations: unperturbed and condition KO
- We describe the distribution across unperturbed cells in gene space as P0
- We describe the distribution across perturbed cells in gene space as Pi
- Further, we are given cell state labels for the unperturbed cells: progenitors, effector, other, cycling, exhausted
- We denote as Q the target cell state proportion vector and Q0 the cell state proportion of the unperturbed cells. And $\deltaQ=Q-Q0$ the desired effect of the perturbation.
- We compute the KNN graph across all cells, that is we include both unperturbed and perturbed cells. For the number of neighbours, we suggest to set K large enough such that neighbourhoods purely consistent of unperturbed cells occur <20% 
- The former proposed number of neighbours K is dependent on the total number of perturbed cells and for good results we advice to exclude condition KO that have fewer than 50 cells in total
- For the statistic, we propose to use neighbourhood enrichment analysis to identify the strength and identity of the perturbation
- More detailed, that means to check in each identified neighbourhood on the KNN graph if the perturbation is positively or negatively enriched
- The enrichment analysis is performed through GLM fits and implemented in the MILO package 
- For a given KO we can count the number of neighbourhoods that have a perturbation signal based on the spatial FDR. We propose a threshold of 0.15.
- Further, we can assign each neighbourhood to a cell state either based on the index cell of the neighbourhood or a majority vote using all unperturbed cells that are included 
- Counting all significant neighbours, we can assign them to one of the cell states.
- This means, we obtain a 5-dimensional vector N that counts the number of neighbourhoods enriched for the perturbation in each one of the cell states, i.e., $N = (P_i, E_i, C_i, Ex_i, O_i)$ being $P_i$ the number of neighbourhoods enriched in the progenitor cell state and $E_i, C_i, Ex_i, O_i$ the equivalent for the other cell states.
Based on the specific number of neighbourhoods identified, it can be advisable to normalise for imbalances. These may occur when one obtains 500 neighbourhoods for progenitors but only 100 for effector states.
- The resulting proportions vector N is our statistic for a given perturbation KO
- To score the perturbation, we further take deltaQ into account. 
- To align our N vector with deltaQ, we count negatively and positively enriched neighbours with -1 and +1, respectively. This is a slight adjustment to the purely positive statistic N that only includes counts. We denote it by D.
- To score the condition KO, we take the L1 distance of the normalised D vector to the deltaQ vector. We denote this as our score S 
- The lower the distance the better the perturbation. If the score is supposed to increase with quality, one could consider to take the inverse 1/S.
