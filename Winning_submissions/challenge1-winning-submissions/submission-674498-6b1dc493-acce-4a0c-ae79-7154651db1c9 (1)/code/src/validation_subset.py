import numpy as np


def min_max_scale(x):
    x = np.array(x, 'float32')
    x = (x - x.min()) / (x.max() - x.min())
    return x


def validation_subset(conditions, labels, size=16):
    subset = []
    # refine based on number of samples of the perturbation
    v = np.array(conditions)
    v, c = np.unique(v, return_counts=True)
    q25, q90 = np.quantile(c, 0.25), np.quantile(c, 0.90)
    for gene, count in zip(v, c):
        if (q25 < count) and (count < q90):
            subset.append(gene)

    # refine based on distance from the mean proportions
    dist_nll, dist_mae = [], []
    mean = np.array((0.0675, 0.2097, 0.3134, 0.3921, 0.0173), 'float32')
    for p in subset:
        y = labels[conditions == p].sum(0)
        y = y / y.sum()
        d = np.abs(y - mean).sum()  # l1
        dist_mae.append(d)
        d = -(mean * np.log(y + 1e-10)).sum()  # neg loglik
        dist_nll.append(d)

    dist_nll, dist_mae = min_max_scale(dist_nll), min_max_scale(dist_mae)
    dist = dist_nll + dist_mae

    # select the 'size' farthest number of genes
    idx = np.argsort(dist)[-size:]
    subset = [subset[i] for i in idx]
    subset = [s for s in subset if s != "Zfp292"]
    return subset

