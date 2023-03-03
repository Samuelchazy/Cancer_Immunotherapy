# Deployment Instructions:

## 1. create environment using conda and activate created environment
$ conda create --name tc_cancer python==3.8

$ conda activate tc_cancer


## 2. install scanpy, scvi-tools, and needed packages
$ conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch -c conda-forge  
$ conda install jax jaxlib cuda-nvcc -c conda-forge -c nvidia  
$ conda install scvi-tools -c conda-forge  

$ conda install -c conda-forge scanpy python-igraph leidenalg  
$ conda install -c conda-forge numpy scipy pandas scikit-learn  
$ conda install -c conda-forge scikit-misc


## 3. execute 'run.py'
$ python run.py

## NOTE: pre-trained model and fresh model training
'run.py' will run by loading the pre-trained model. 
If you want to train a model fresh to repeat the result, uncomment 3-lines in the model training block in 'run.py'.  Please, refer to the 'run.py' for details.