
# coding: utf-8

import sys

import pandas as pd



df = pd.read_csv(sys.argv[1])
df.head()


samples = [s[:-2] for s in df.columns[2::4]]


with open(sys.argv[2], "w") as f:
    f.write(",median\n")
    for s in samples:
        cov = df[[s+'-'+N for N in ['A', 'C', 'G', 'T']]].sum(axis=1)
        f.write(s + ',' + str(cov.median()) + '\n')
        



