import argparse

import matplotlib
matplotlib.use("AGG")

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

import seaborn as sns
import pandas as pd

from sklearn.decomposition import PCA

sns.set_context('poster')
sns.set_style('white')

pd.options.mode.chained_assignment = None  # default='warn'

import hdbscan
from collections import Counter
from collections import defaultdict
from numpy import random

#_____________________________


def normalize(x, r):
    M = np.divide(x, r)
    M_norm = np.full_like(M, 0)
    for i in range(np.shape(M)[0]):
        rev = 1 - M[i, :]
        if np.dot(M[i, :], M[i, :]) > np.dot(rev, rev):
            M_norm[i, :] = rev
        else:
            M_norm[i, :] = M[i, :]
    return M_norm


def draw_PCA(f_pca, black_points=None):
    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel("PC 1")
    plt.ylabel("PC 2")

    plt.scatter(f_pca[:, 0], f_pca[:, 1], s=10, linewidth=0, alpha=0.2);
    if black_points is not None:
        plt.scatter(f_pca[black_points, 0], f_pca[black_points, 1], s=10, linewidth=0, c="black", alpha=1);
        plt.title("%s/%s points" % (np.sum(black_points), len(f_pca)))
    else:
        plt.title("%s points" % len(f_pca))
    plt.savefig("logs/all_SNPs.png")


def filter_by_coverage(cur_r, bad_percent, bad_samples):
    def filter_row(row):
        num_of_samples = len(row)
        valid = np.sum(np.array(([(min_coverage < row) & (row < max_coverage)])))
        return num_of_samples - valid <= bad_samples

    min_coverage = np.percentile(cur_r, bad_percent, axis=0)
    max_coverage = np.percentile(cur_r, 100-bad_percent, axis=0)
    good_coverage = np.array([filter_row(row) for row in cur_r])
    return good_coverage


def draw_legend(class_colours, classes, right=False):
    recs = []
    for i in range(0, len(classes)):
        recs.append(mpatches.Rectangle((0,0), 1, 1, fc=class_colours[i]))
    if right:
        plt.legend(recs, classes, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        plt.legend(recs, classes)

#____________________________________________________________________________________________
def sizes_to_probs(s):
    s_log = np.log(s)
    return s_log / s_log.sum()


def subsample_from_cluters(Npoints, labels):
    groupby_labels = defaultdict(list)
    for i in range(len(labels)):
        groupby_labels[labels[i]].append(i)
    groupby_labels = dict(groupby_labels)

    cluster_size = np.array([len(groupby_labels[i]) for i in groupby_labels.keys()], dtype=int)
    print("\nCluster sizes:")
    print(list(cluster_size), "\n")
    cluster_size = Npoints * sizes_to_probs(cluster_size)
    print("Subsamplings sizes:")
    print([int(round(s)) for s in cluster_size], "\n")

    subsamples = {}
    for i in groupby_labels.keys():
        subsample_size = int(round(cluster_size[i]))
        subsample_index = random.choice(groupby_labels[i], size=subsample_size, replace=True).tolist()
        subsamples[i] = subsample_index

    ind = [val for sublist in list(subsamples.values()) for val in sublist]
    random.shuffle(ind)
    return ind


def clusterization(f, x, r, x_new, r_new, how_many_remain, pca=True, num_of_comp=2):
    if pca:
        f_pca = PCA(n_components=num_of_comp).fit(f).transform(f)
        cur_f = f_pca
    else:
        cur_f = f
        f_pca = PCA(n_components=2).fit(f).transform(f)

    # N = (nt) (len(f) * 0.005)
    # print(N)
    N = 100

    clusterer = hdbscan.HDBSCAN(min_cluster_size=N, min_samples=1).fit(cur_f)

    plt.figure(figsize=(15, 8))
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel("PC 1")
    plt.ylabel("PC 2")
    if pca:
        plt.title("Clustering %s primary components" % num_of_comp)
    else:
        plt.title("Clustering initial frequencies")

    color_palette = sns.color_palette("Set2", 20)
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in clusterer.labels_]
    cluster_member_colors = [sns.desaturate(x, p) for x, p in
                             zip(cluster_colors, clusterer.probabilities_)]
    plt.scatter(f_pca[:, 0], f_pca[:, 1], s=40, linewidth=0, c=cluster_member_colors, alpha=0.3);

    not_outliers = np.where(clusterer.labels_ != -1)[0]
    f_no_outliers = f[not_outliers, :]
    x_no_outliers = x[not_outliers, :]
    r_no_outliers = r[not_outliers, :]
    labels_no_outliers = clusterer.labels_[not_outliers]

    subs_ind = subsample_from_cluters(how_many_remain, labels_no_outliers)

    plt.scatter(f_pca[:, 0], f_pca[:, 1], s=10, linewidth=0, c=cluster_colors, alpha=0.3)
    plt.scatter(f_pca[subs_ind, 0], f_pca[subs_ind, 1], s=10, linewidth=0, c="black", alpha=1)
    plt.savefig("logs/clusters.png")

    np.savetxt(r_new, r_no_outliers[subs_ind, :], fmt='%s')
    np.savetxt(x_new, x_no_outliers[subs_ind, :], fmt='%s')

    sizes_of_classes = Counter(clusterer.labels_)
    print(sizes_of_classes.get(-1, 0), "outliers\n")
    labels = [str(x) + ' - ' + str(sizes_of_classes[x]) for x in range(max(clusterer.labels_) + 1)]
    draw_legend(color_palette, labels, right=True)

    print("Medians in clusters:")
    for i in range(max(clusterer.labels_) + 1):
        f_with_labels = f.copy()
        f_with_labels = np.hstack([f_with_labels, clusterer.labels_.reshape(len(f_with_labels), 1)])
        col = f_with_labels[:, -1]
        idx = (col == i)
        print(i, np.round(np.median(f_with_labels[idx, :-1], axis=0), 2))

    plt.savefig("logs/clusterization.png")

#________________________________________________________________________________________________________


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("r_old", type=str,
                        help="Initial R matrix")
    parser.add_argument("x_old", type=str,
                        help="Initial X matrix")
    parser.add_argument("r_new", type=str,
                        help="Filtered R matrix")
    parser.add_argument("x_new", type=str,
                        help="Filtered X matrix")

    parser.add_argument("-k", "--num_of_filtered_snps", type=int, default=200,
                        help="How many SNPs to ramain")

    parser.add_argument("-p", "--num_of_pure_samples", type=int, default=0,
                        help="Number of pure samples in data")

    parser.add_argument("-norm", "--normalize", action="store_true",
                        help="Normalize frequencies is this flag specified")

    parser.add_argument("-cp", "--coverage_percentile", type=int, default=15,
                        help="What percentiles to choose for coverage filtering")

    parser.add_argument("-bs", "--num_of_bad_samples", type=int, default=2,
                        help="How many samples to ignore during coverage filtering")

    args = parser.parse_args()

    r = np.genfromtxt(args.r_old, dtype=int, delimiter=' ')
    x = np.genfromtxt(args.x_old, dtype=int, delimiter=' ')

    print("%s sites" % len(r))

    r = np.delete(r, range(args.num_of_pure_samples), axis=1)
    x = np.delete(x, range(args.num_of_pure_samples), axis=1)

    Ncut = 6
    print("\nDelete zero and almost zero profiles:")
    good_ind = [i for i in range(np.shape(x)[0])
                if not ((np.abs(r[i, :] - x[i, :]) <= Ncut).all() or (x[i, :] <= Ncut).all())]
    print(len(good_ind), "remained")

    x = x[good_ind, :]
    r = r[good_ind, :]

    if args.normalize:
        f = normalize(x, r)
    else:
        f = np.divide(x, r)

    f_pca = PCA(n_components=2).fit(f).transform(f)

    good_coverage = filter_by_coverage(r, args.coverage_percentile, args.num_of_bad_samples)
    draw_PCA(f_pca, good_coverage)

    f = f[good_coverage, :]

    clusterization(f, x, r, args.x_new, args.r_new, args.num_of_filtered_snps, pca=True, num_of_comp=5)


if __name__ == '__main__':
    main()


