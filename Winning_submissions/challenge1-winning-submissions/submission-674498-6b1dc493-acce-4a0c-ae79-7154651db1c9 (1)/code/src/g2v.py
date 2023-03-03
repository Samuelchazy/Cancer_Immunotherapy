import numpy as np
import wget
import os


G2V_URL = "https://github.com/jingcheng-du/Gene2vec/raw/master/pre_trained_emb/gene2vec_dim_200_iter_9_w2v.txt"
EMB_PATH = "./data/embeddings/gene2vec_dim_200_iter_9_w2v.txt"


def load_embeddings():
    if not os.path.isfile(EMB_PATH):
        wget.download(G2V_URL, EMB_PATH)

    with open(EMB_PATH, 'r') as f:
        lines = f.readlines()[1:]

    g2v_embeddings = {}
    for x in lines:
        x = x.split(' ')
        g = x[0].upper()
        g2v_embeddings[g] = np.asarray(x[1:-1], 'float32')

    return g2v_embeddings

