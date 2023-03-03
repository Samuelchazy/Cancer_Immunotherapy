#####
# The sc_training.h5ad file needs to be in the same directory as the run.py script
# The resulting files will be produced in the same directory as the script
# Run the following commands to set up the appropriate environment
# All needed files to run the script are located in the Google Drive link below:
# https://drive.google.com/drive/folders/1x2g0XhNBWU6saf77Z4xgeobaA3xEXoUI?usp=share_link
#####

conda create -n topcoder

conda activate topcoder

conda install -c fastai fastai

conda install -c conda-forge scanpy python-igraph leidenalg

python3 run.py



