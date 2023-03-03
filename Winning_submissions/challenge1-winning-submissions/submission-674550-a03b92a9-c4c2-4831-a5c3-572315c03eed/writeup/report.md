## 1-page writeup
My solution is based on the assumption that there would exist a subgroup of 66 KO genes that have a significant correlation with the target gene in the latent space, and the cell state distribution of the selected sub-group genes could proportionaly reflect those of the target gene.  
So, the steps that I had taken to generate the submitted solutions are:  

1. Construct a prior-based naive model as a baseline:  
    The validation and test genes' KO cell state distritution have been calulated based on a naive assumption that when a gene has 'zero' expression value in 'unperturbated' cell, it has been considered as an instance of knock-out of that gene as a prior.  
So, once collecting all such cases (i.e. all the cells that show zero expression value of the corresponding gene in unperturbed condition) as the hypothetical knock-out cases, the cell state distribution has been calculated in those cells coming from 'unperturbed' condition.

2. Choose a model to represent latent space of the cell populations of all(or some) conditions:  
    After literature survey, [this paper](https://pubmed.ncbi.nlm.nih.gov/32176273/) by Svensson et al.[1] drew my attention. I followd the lead of theirs and applied  linearly decoded variational autoencoder(LDVAE) to the given dataset to obtain Z-representation using [scvi-tools](https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html)[2].

3. Select a creterion to measure the correlation of the 66KO genes with the target gene in the latent space:  
    In the linear decodable latent space obtained From step2, weights of latent variables determine the variation of the gene in that space, so the correlations of weight vectors among genes could suggest that the genes with high correlation could possibly co-express. Hence, Pearson correlation coefficient had been calculated between the latent varialbe weights of target genes and those of 66 KO genes, and the genes with pearson-r > 0.7 and p-value <0.05 had been selected as the 'highly correlated' subgroup of 66 KO genes. see 'run.py' for the details. 


4. Calculate the cell state distribution by counting cell number of each state of the selected subgroup of KO genes:  
    see 'run.py' for the details


5. Build the solution by mixing the baseline solution from step1 and the gene correlation based solution from step4:  
    Here, the baseline solution functions as both a regulator and a fallback choice where there has not been found any signficantly correlated KO gene, and amount of mixing can be controlled by 'alpha' (=0: baseline solution, 1: correlation based solution). see 'run.py' for the details.

Following these steps, in a late stage of the competition, I have tested my assumption by setting up my local validation condtion, which was to subset training dataset with arbitrarily chosen validation genes = ['Tox2', 'Arid4b', 'Tcf7', 'Il12rb2'], and observed that l1-loss have improved from 1.1801 (nailve baseline solution) to 0.8840 (gene correlation based solution)

REFERENCES  
[1] Svensson,V. et al. (2020) Interpretable factor models of single-cell RNA-seq
via variational autoencoders. Bioinformatics, 36, 3418–3421  

[2] Linearly-decoded Variational Auto-encoder (LDVAE), https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html, accessed on 01/27/2023 