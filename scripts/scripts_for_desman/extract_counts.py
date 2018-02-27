
# coding: utf-8

import pandas as pd
import os, re, glob
import sys


dir_with_counts = sys.argv[1]
output_file = sys.argv[2]


S = 0
counts = []
for name in glob.glob(dir_with_counts+"/*.cnt"):
    basename = os.path.basename(name)    
    m = re.search("(.*)\.cnt",basename)
    if m:
        sample_name = m.groups()[0] 
        S = S + 1
        counts.append((name, sample_name))

all_dfs = []
col_names = ["Contig", "Position"]
for file, sample_name in counts:
    df = pd.read_csv(file, sep='\t', 
                     names=["contig", "pos", "ref", "depth", "all", "A", "C", "G", "T", "N"])
    df.drop(["ref", "depth", "all", "N"], axis=1, inplace=True)
    
    for base in ["A", "C", "G", "T"]:
        df[base] = df[base].str.split(pat=":", expand=True)[1]
        col_names.append(sample_name+"-"+base)
    
    all_dfs.append(df)


cur_merge = all_dfs[0]

for i in range(1, S):
    cur_merge = pd.merge(cur_merge, all_dfs[i],  how='inner', 
                         left_on=["contig", "pos"], right_on = ["contig", "pos"])

cur_merge.columns = col_names

cur_merge.to_csv(output_file, index=False)

