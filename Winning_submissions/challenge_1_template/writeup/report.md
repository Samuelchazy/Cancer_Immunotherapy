Challenge submission:
Participant: huyennguyen29 (Thu Huyen Nguyen)

The prediction model include 2 parts:

1. The pretrained gene embedding (vector of 100 dimensions) is from scBERT paper (https://github.com/TencentAILabHealthcare/scBERT). The code to extract gene embedding vector for all genes can be found in the gen2vec.ipynb notebook. The gen2vec weights are stored in ./data folder (it is the original weights taken from its github, no finetuning is performed). However, the pretrained model was used human genome, therefore, an example of scRNA-seq data is provided (./data/panglao_10000.h5ad) to match the genenames provided by the challenge.

2. Perform a GridSearch for KNeighborsRegressor and RandomForestRegressor to map gene embedding vector to 5-state vectors. A StandardScaler step on all gene embedding was performed before training the regression models. 