import matplotlib
matplotlib.use("AGG")

import numpy as np
import pandas as pd
import seaborn as sns
import sys

import matplotlib.pyplot as plt

import hdbscan
from collections import Counter
from sklearn.decomposition import PCA

from numpy import random
from collections import defaultdict

sns.set_context('poster')
sns.set_style('white')
sns.set_color_codes()
plot_kwds = {'alpha': 0.2, 's': 10, 'linewidths': 0}


def normalize(M):
    M_norm = np.full_like(M, 0)
    for i in range(np.shape(M)[0]):
        rev = 1 - M[i, :]
        if np.dot(M[i, :], M[i, :]) > np.dot(rev, rev):
            M_norm[i, :] = rev
        else:
            M_norm[i, :] = M[i, :]

    return M_norm


def sizes_to_probs(s):
    s_log = np.log(s)
    return s_log / s_log.sum()


def subsample_from_cluters(Npoints, data, labels):
    groupby_labels = defaultdict(list)
    for i in range(np.shape(data)[0]):
        groupby_labels[labels[i]].append(i)
    groupby_labels = dict(groupby_labels)

    cluster_size = np.array([len(groupby_labels[i]) for i in groupby_labels.keys()], dtype=int)
    cluster_size = Npoints * sizes_to_probs(cluster_size)
    print("\nSubsamplings sizes:")
    print(cluster_size)

    subsamples = {}
    for i in groupby_labels.keys():
        subsample_size = int(round(cluster_size[i]))
        subsample_index = random.choice(groupby_labels[i], size=subsample_size, replace=True).tolist()
        subsamples[i] = subsample_index

    ind = [val for sublist in list(subsamples.values()) for val in sublist]
    random.shuffle(ind)
    return ind


def filter_matrices(r_name, x_name, Npoints=100, Ncut=5):
    r = np.genfromtxt(r_name, dtype=int, delimiter=' ')
    x = np.genfromtxt(x_name, dtype=int, delimiter=' ')

    print("\nBefore filtration:")
    print(np.shape(x)[0])

    print("\nDelete zero and almost zero profiles:")
    good_ind = [i for i in range(np.shape(x)[0]) if
                not ((np.abs(r[i, :] - x[i, :]) <= Ncut).all() or (x[i, :] <= Ncut).all())]
    print(len(good_ind))

    x = x[good_ind, :]
    r = r[good_ind, :]
    f = normalize(np.divide(x, r))

    N_min_cluster = (int)(np.shape(f)[0] * 0.01)
    clusterer = hdbscan.HDBSCAN(min_cluster_size=N_min_cluster).fit(f)

    print("\nDelete outliers:")
    clusters = Counter(clusterer.labels_)
    print(np.shape(f)[0] - clusters[-1])
    print("\nCLusters (-1 <-> outliers):")
    print(clusters)

    not_outliers = np.where(clusterer.labels_ != -1)[0]
    f_no_outliers = f[not_outliers, :]
    x_no_outliers = x[not_outliers, :]
    r_no_outliers = r[not_outliers, :]
    labels_no_outliers = clusterer.labels_[not_outliers]
    f_pca = PCA(n_components=2).fit(f_no_outliers).transform(f_no_outliers)

    color_palette = sns.color_palette("Set2", 20)
    cluster_colors = np.array([color_palette[x] for x in labels_no_outliers])

    subs_ind = subsample_from_cluters(Npoints, f_no_outliers, labels_no_outliers)

    plt.scatter(f_pca[:, 0], f_pca[:, 1], s=10, linewidth=0, c=cluster_colors, alpha=0.3)
    plt.scatter(f_pca[subs_ind, 0], f_pca[subs_ind, 1], s=10, linewidth=0, c="black", alpha=1)
    plt.savefig("logs/clusters.png")

    np.savetxt(sys.argv[3], r_no_outliers[subs_ind, :], fmt='%s')
    np.savetxt(sys.argv[4], x_no_outliers[subs_ind, :], fmt='%s')


filter_matrices(sys.argv[1], sys.argv[2], int(sys.argv[5]))
