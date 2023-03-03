import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import mahalanobis

def score(P0, Q = (0.6, 0.25, 0.05, 0.1, 0.00), Pi):
    inv_cov = np.linalg.inv(np.cov(P0.T))
    mahal_dist = np.sqrt(mahalanobis(Pi, P0, inv_cov) ** 2)
    cos_sim = cosine_similarity(Q.reshape(1, -1), Pi.reshape(1, -1))[0][0]
    return (1 - mahal_dist) * cos_sim
