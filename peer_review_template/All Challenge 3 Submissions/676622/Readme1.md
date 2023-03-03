## CHALLANGE 1 is in "dataset.py"

The dataset is formated anndata in the file sc_training.h5ad

total counts per cell followed by log(1+x)-transformation

the raw gene expression values before normalization in .layers[‘rawcounts’].

The condition information ('Unperturbed' or the target gene name if perturbed) is provided in .obs[‘condition’]
The cell state of each single cell is stored in .obs[‘state’]

## Requirements ##
Predict the cell-state proportions of the held-out genes: For each of the 7 held-out knockouts (targeting genes 'Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3')
Cell state proportions should add to be a+b+c+d+e=1.

## Result ##
validation_output.csv: This file should contain the 3×5 matrix in csv format, as well as the name of the gene in the first column.
test_output.csv: This file should contain the 4×5 matrix in csv format.  the 5-demensional vector prediction for the 4 held-out test set genes.