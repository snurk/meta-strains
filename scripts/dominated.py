
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import os
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
import sys


def p2f(x):
    # percents to float
    if x == '-':
        return None
    else:
        return float(x.strip('%'))/100


def find_margin(VAFs, sample_name=None, within_margin=0.05, eps=0.01):
    
    print("Num of SNVs:", len(VAFs))
        
    margin = eps
    while ((0.5 - margin < VAFs) & (VAFs < 0.5 + margin)).sum() < within_margin * len(VAFs) and 0.5 > margin:
        margin += eps

    if margin > 0.5:
        margin = 0.5

    print('margin = %.2f' % margin)

    if 0.5 + margin >= 0.7:
        color = 'red'
        print('DOMINATED')
        res = True
    else:
        color = 'blue'
        print('NOT dominated')
        res = False
        
    print("%.2f - %.2f" % (0.5-margin, 0.5+margin))

    plt.hist(VAFs, 50, alpha=0.5, color=color);
    plt.xlim((0, 1))
    plt.xlabel('SNV frequencies')
    plt.title(sample_name)
    plt.savefig("hists/%s.png" % sample_name)
    
    return(res)


def filter_by_coverage(depth, vafs, bad_percentile=0.3, good_samples_percent=0.8):
    q1 = depth.quantile(bad_percentile)
    q2 = depth.quantile(1 - bad_percentile)
    
    cur_n_samples = depth.shape[1]
    necessary_amount = int(cur_n_samples * good_samples_percent)
    
    ind = ((depth > q1) & (depth < q2)).sum(axis=1) >= necessary_amount
    
    return ind



def main():
    df = pd.read_csv(sys.argv[1], sep='\t')

    df_samples = df.iloc[:,-1].str.split(pat=" ", expand=True)
    n_samples = df_samples.shape[1]

    df_samples_cov = df_samples.apply(lambda x : x.str.split(pat=":", expand=True)[1]).astype("int64")
    df_samples_VAF = df_samples.apply(lambda x : x.str.split(pat=":", expand=True)[4]).applymap(p2f)

    sample_names = sys.argv[2].split(',') # ["sample"+str(i) for i in range(1, 12)]

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
    genotypes = genotypes[(genotypes.sum(axis=1) > 0) & (genotypes.sum(axis=1) < len(dominated_samples))] #remove non-informative sites

    g = sns.clustermap(genotypes.T, 
                       xticklabels = False, 
                       yticklabels=[sample_names[i] for i in dominated_samples])
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.cax.set_visible(False)
    plt.suptitle('%i SNVs' % len(genotypes))
    plt.savefig("dominated_genotypes.png")


if __name__ == "__main__":
    main()


