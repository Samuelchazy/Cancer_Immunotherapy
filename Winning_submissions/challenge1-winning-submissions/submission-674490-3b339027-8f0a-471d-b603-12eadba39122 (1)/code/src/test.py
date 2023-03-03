import argparse
import pandas as pd

from joblib import load
from utils import load_features


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--jsonl-path', required=True,
                        help='path to the jsonl file used to build features.')
    parser.add_argument('--model-path', required=True,
                        help='path to the model.')
    parser.add_argument('--genes', required=True, nargs='+',
                        help='list of genes to predict.')
    parser.add_argument('--csv-path', required=True,
                        help='path where the solution will be saved.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parser()
    if len(args.genes) == 0:
        raise Exception('No gene supplied for prediction.')

    # Load the model and the features
    with open(args.model_path, 'rb') as f:
        model = load(f)
    features = load_features(args.jsonl_path)

    # Make predictions
    df = pd.DataFrame(
        index=args.genes,
        columns=['a_i', 'b_i', 'c_i', 'd_i', 'e_i']
    )
    df.index.name = 'gene'
    df.loc[args.genes] = model.predict(features.loc[args.genes].to_numpy())
    df = df.div(df.sum(axis=1), axis=0)  # probabilities sum to 1

    # Save
    df.to_csv(args.csv_path, index=True)
    print(f'Predictions saved into {args.csv_path}.')
