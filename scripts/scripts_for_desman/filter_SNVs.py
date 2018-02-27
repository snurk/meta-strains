
# coding: utf-8

import pandas as pd
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]
num_SNVs = int(sys.argv[3])


df = pd.read_csv(input_file, sep='\t')

df_samples = df.iloc[:,-1].str.split(pat=" ", expand=True)

df_samples_cov = df_samples.apply(lambda x : x.str.split(pat=":", expand=True)[1]).astype("float64")

cov_median = df_samples_cov.median()

selected_indexes = ((df_samples_cov.subtract(cov_median, axis=1)).pow(2)).sum(axis=1).sort_values().index[:num_SNVs]

selected_sites = df.loc[selected_indexes, ["Chrom", "Position", "Position"]]

selected_sites.to_csv(output_file, sep='\t', header=False, index=False)
