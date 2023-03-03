import numpy as np


def sample_dist(predict_func, perturbation_emb, n_samples):
    p = np.array([perturbation_emb for _ in range(n_samples)])
    z = np.random.normal(scale=1., size=(n_samples, 64)).astype('float32')

    dist = predict_func(p, z)[0]
    dist = dist.numpy().sum(0)

    dist = dist / dist.sum()
    # dist[dist < 0.02] = 0
    dist = dist.round(3)
    dist = dist / dist.sum()

    return dist


def evaluate(predict_func, perturbations_test, condition_test, label_test,
             g2v_embeddings):
    scores = []

    for p in perturbations_test:
        y = label_test[condition_test == p]
        y = y.sum(0)
        true_proportions = y / y.sum()

        q = np.zeros((200,), 'float32')
        if p.upper() in g2v_embeddings:
            q = g2v_embeddings[p.upper()]
        pred_proportions = sample_dist(predict_func, q, 1000)

        s = np.abs(pred_proportions - true_proportions).sum()
        scores.append(s)

    return np.mean(scores)

