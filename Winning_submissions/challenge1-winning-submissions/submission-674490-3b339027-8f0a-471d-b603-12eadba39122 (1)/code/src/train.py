import argparse
import numpy as np
import pandas as pd
import scanpy as sc

from sklearn import ensemble
from joblib import dump
from utils import load_features


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--jsonl-path', required=True,
                        help='path to the jsonl file used to build features.')
    parser.add_argument('--h5-path', required=True,
                        help='path to the h5 file containing the targets.')
    parser.add_argument('--model-path', required=True,
                        help='path where the model will be saved.')
    parser.add_argument('--min-cells', required=True, type=int, default=25,
                        help='minimum number of cells for a condition to be considered.')
    return parser.parse_args()


def load_targets(path_to_sc, min_cells_per_cond):
    """Generate the targets

    Args:
        path_to_sc: path to the h5 file
        min_cells_per_cond: only conditions with at least this number of cells are considered
    """

    adata = sc.read_h5ad(path_to_sc)
    n_cells_per_cond = adata.obs.groupby('condition').size()
    excluded_conds = n_cells_per_cond[n_cells_per_cond <
                                      min_cells_per_cond].index.tolist()
    if excluded_conds:
        print(f'Conditions without enough cells: {excluded_conds}')

    targets = (
        adata.obs
        .groupby('condition')['state']
        .value_counts(normalize=True)
        .unstack()
        .drop(excluded_conds + ['Unperturbed'])
        .loc[:, ['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other']]
    )

    return targets


if __name__ == '__main__':
    args = parser()

    # Load the data
    features = load_features(args.jsonl_path)
    targets = load_targets(args.h5_path, args.min_cells)

    # Define the model
    model = ensemble.ExtraTreesRegressor(
        n_estimators=100,
        criterion='absolute_error',
        max_leaf_nodes=2,
        n_jobs=-1
    )

    # Train
    model = model.fit(
        features.loc[targets.index].to_numpy(),
        targets.to_numpy()
    )

    # Save
    with open(args.model_path, 'wb') as f:
        dump(model, f)
    print(f'Model saved into {args.model_path}')
