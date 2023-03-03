## Installation

The required packages can be installed with the following command.
```
pip install -r requirements.txt
```
It is advised to first create a virtual environment.

## Training
The `code/train.sh` script shows how to call the training program. You first need to download the h5 file containing the targets. The features are supplied in the `data` folder). Also, conditions which have less than the supplied minimum number of cells will be ignored from training.
```sh
python src/train.py \
  --jsonl-path data/data_report.jsonl \
  --h5-path data/sc_training.h5ad \
  --model-path checkpoint/model.joblib \
  --min-cells 25
```

## Predictions
The `code/test.sh` script shows how to call the training program. You need to supply a list of genes and indicate where the model is saved.
```sh
python src/test.py \
  --jsonl-path data/data_report.jsonl \
  --model-path checkpoint/model.joblib \
  --genes "Ets1" "Fosb" "Mafk" "Stat3" \
  --csv-path ../solution/test_output.csv
```