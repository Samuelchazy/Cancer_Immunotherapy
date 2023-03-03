import jsonlines
import numpy as np
import pandas as pd


def select_transcript(transcripts):
    """Select a single transcript among  a list"""
    candidates = [t for t in transcripts
                  if t['type'] == 'PROTEIN_CODING'
                  and 'exons' in t]
    n_candidates = len(candidates)
    if n_candidates == 0:
        raise Exception('No transcript found')

    if n_candidates == 1:
        # Only 1 transcript => no selection needed
        return candidates[0]

    # Select transcript by name
    for candidate in candidates:
        if candidate.get('name', '').lower().strip() == 'transcript variant 1':
            return candidate

    # If everything failed, return the first transcript
    return candidates[0]


def load_features(path_to_report):
    """Generate the features

    Args:
        path_to_report: path to the jsonl file
    """
    features = dict()
    with jsonlines.open(path_to_report) as reader:
        for obj in reader:
            sym = obj['symbol']
            transcript = select_transcript(obj['transcripts'])
            exons_len = [int(e['end']) - int(e['begin']) +
                         1 for e in transcript['exons']['range']]
            sym_feats = dict()
            sym_feats['orientation'] = - \
                1 if obj['orientation'] == 'minus' else 1
            sym_feats['exons_count'] = len(exons_len)
            sym_feats['exons_len_mean'] = np.mean(exons_len)
            sym_feats['exons_len_min'] = np.min(exons_len)
            sym_feats['exons_len_max'] = np.max(exons_len)
            sym_feats['exons_len_median'] = np.median(exons_len)
            features[sym] = sym_feats
    return pd.DataFrame.from_dict(features, orient='index')
