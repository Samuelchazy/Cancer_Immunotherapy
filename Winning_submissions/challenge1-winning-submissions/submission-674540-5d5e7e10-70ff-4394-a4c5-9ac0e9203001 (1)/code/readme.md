#### Deployment Instructions

Our code and data are shared in https://www.dropbox.com/sh/ugid5xma7x7sxdk/AAAG1OgEb4AZEc_VZD90kB6la?dl=0.

##### Code files
* The code for generating the prediction vectors is in `code/human_to_mouse.ipynb`. Running the whole notebook would produce our submitted output. To make prediction for new genes, modify `knockouts_test` variable.

##### Data files
* `raw_data/sc_training.h5ad` is the given data file.
* `processed_data/mouse_adata_pca.h5ad` is the post-PCA result of `raw_data/sc_training.h5ad`. That is, we run PCA with 50 dimension, and replace adata.X with the post-PCA results.
* `processed_data/human_adata_pca_ver2.h5ad` is the post-PCA result of public data set `K562_gwps_normalized_singlecell_01.h5ad` in `https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387`.
* `processed_data/mouse_to_human.txt` is a dictionary for translating mouse gene names to human gene names. Most of the time the translation is just taking upper case letters.
