import pandas as pd
import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt
import os
import seaborn as sns
import warnings

import scipy
import scipy.cluster

warnings.filterwarnings('ignore')
import sys


def p2f(x):
    # percents to float
    if x == '-':
        return None
    else:
        return float(x.strip('%')) / 100


def find_margin(VAFs, sample_name=None, within_margin=0.05, eps=0.01):
    print("Num of SNVs:", len(VAFs))

    margin = eps
    while ((0.5 - margin < VAFs) & (VAFs < 0.5 + margin)).sum() < within_margin * len(VAFs) and 0.5 > margin:
        margin += eps

    if margin > 0.5:
        margin = 0.5

    if len(VAFs) < 500:
        print("Very few SNVs")
    else:
        print("%.2f - %.2f" % (0.5 - margin, 0.5 + margin))
        print('margin = %.2f' % margin)

    if 0.5 + margin >= 0.7 or len(VAFs) < 500:
        color = 'red'
        print('DOMINATED')
        res = True
    else:
        color = 'blue'
        print('NOT dominated')
        res = False

    plt.hist(VAFs, 50, alpha=0.5, color=color)
    plt.xlim((0, 1))
    plt.xlabel('SNV frequencies')
    plt.title(sample_name)
    plt.savefig("hists/%s.png" % sample_name)

    return res


def filter_by_coverage(depth, vafs, bad_percentile=0.3, good_samples_percent=0.8):
    q1 = depth.quantile(bad_percentile)
    q2 = depth.quantile(1 - bad_percentile)

    cur_n_samples = depth.shape[1]
    necessary_amount = int(cur_n_samples * good_samples_percent)

    ind = ((depth > q1) & (depth < q2)).sum(axis=1) >= necessary_amount

    return ind


def draw_heatmap(genotypes, sample_names, clusters, dominated_samples):
    colors = ['yellow', 'purple', 'orange', '#96cde6', 'red', '#c0bd7f', 'green', '#5fa641', '#d485b2',
              '#4277b6', '#df8461', '#463397', '#e1a11a', '#91218c', '#e8e948', '#7e1510',
              '#92ae31', '#6f340d', '#d32b1e', '#2b3514']

    last_color = 0
    row_colors = [-1] * len(sample_names)
    for clustered_samples in clusters.values():
        for sample in clustered_samples:
            row_colors[sample_names.index(sample)] = colors[last_color]
        last_color += 1
    row_colors = [row_colors[i] for i in dominated_samples]

    g = sns.clustermap(genotypes,
                       xticklabels=False,
                       yticklabels=[sample_names[i] for i in dominated_samples],
                       cmap="plasma",
                       row_colors=row_colors)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.cax.set_visible(False)
    plt.suptitle('%i SNVs' % genotypes.shape[1])
    #plt.interactive(False)
    #plt.show(block=True)
    plt.savefig("dominated_genotypes.png")


def main():
    df = pd.read_csv(sys.argv[1], sep='\t')

    df_samples = df.iloc[:, -1].str.split(pat=" ", expand=True)
    n_samples = df_samples.shape[1]

    df_samples_cov = df_samples.apply(lambda x: x.str.split(pat=":", expand=True)[1]).astype("int64")
    df_samples_VAF = df_samples.apply(lambda x: x.str.split(pat=":", expand=True)[4]).applymap(p2f)

    sample_names = sys.argv[2].split(',')

    if not os.path.exists("hists"):
        os.makedirs("hists")

    dominated_samples = []
    for i in range(n_samples):
        print(sample_names[i])

        cur_depth = df_samples_cov.iloc[:, i].copy()
        q1 = cur_depth.quantile(0.3)
        q2 = cur_depth.quantile(0.7)

        print("Coverage median: %i" % cur_depth.median())
        print("Q30:             %i" % q1)
        print("Q70:             %i" % q2)

        cur_VAFs = df_samples_VAF.iloc[:, i]
        selected_VAFs = cur_VAFs[(q1 < cur_depth) & (cur_depth < q2)]
        selected_VAFs = selected_VAFs[(selected_VAFs > 0.02) & (selected_VAFs < 0.98)]

        plt.figure(i)

        if find_margin(selected_VAFs, sample_name=sample_names[i]):
            dominated_samples.append(i)

        print()

    df_dominated_cov = df_samples_cov.iloc[:, dominated_samples]
    df_dominated_VAF = df_samples_VAF.iloc[:, dominated_samples]

    selected_SNVs = filter_by_coverage(df_dominated_cov, df_dominated_VAF)

    # clustering of genotypes in dominated samples
    genotypes = df_dominated_VAF[selected_SNVs] > 0.5
    # genotypes = genotypes[
    #    (genotypes.sum(axis=1) > 0) & (genotypes.sum(axis=1) < len(dominated_samples))]  # remove non-informative sites
    genotypes = genotypes.T

    # percent of non-matching SNVs
    dists = scipy.spatial.distance.pdist(genotypes, 'hamming')

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    Z = scipy.cluster.hierarchy.linkage(dists, 'complete')

    clusters = {}
    for i in range(len(dominated_samples)):
        clusters[i] = {sample_names[dominated_samples[i]]}
    last_cluster = len(dominated_samples) - 1

    for i in range(len(Z)):
        if Z[i, 2] > 0.01:
            break

        u, v = Z[i, 0], Z[i, 1]
        last_cluster += 1
        clusters[last_cluster] = clusters[u] | clusters[v]  # union
        clusters.pop(u)
        clusters.pop(v)

    print("Clustering results:")
    for clustered_samples in clusters.values():
        print(", ".join(clustered_samples))

    draw_heatmap(genotypes, sample_names, clusters, dominated_samples)


if __name__ == "__main__":
    main()
