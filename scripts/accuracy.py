import pandas as pd
import numpy as np

from collections import Counter
from itertools import permutations

import networkx as nx


def print_edges(edges):
    print(','.join([str(e) for e in edges]))


# Считываем граф
# Oставляем только большую компоненту связности
G = nx.DiGraph()
with open("assembly/K127/prelim_graph.gfa") as f:
    for line in f:
        line = line.split()
        line_type = line[0]
        
        # S 238024 ACCAATTAT KC:i:37210
        if line_type == "S":
            v_name = int(line[1])
            v_length = len(line[2])
            G.add_node(v_name, length=v_length)
        
        # L 238322 + 19590 - 55M
        if line_type == "L":
            v1 = int(line[1])
            v2 = int(line[3])
            G.add_edge(v1, v2)


# remain only largest component
new_G = nx.DiGraph()
for g in nx.weakly_connected_component_subgraphs(G):
    #print(g.number_of_nodes())
    if new_G.number_of_nodes() < g.number_of_nodes():
        new_G = g.copy()
G = new_G.copy()


# Табличка с референсами
# Считываем файл ответа, как он есть
df_ref = pd.read_csv("refs/refs_edges.txt", header=None, names=["e"])
df_ref = df_ref["e"].str.split('\t', 1, expand=True)
df_ref.columns = ["e_id", "strains"]
df_ref = df_ref.set_index("e_id")
df_ref.index = df_ref.index.astype("int")


# Оставляем только ребра из большой компоненты:
df_ref = df_ref.loc[list(G.nodes)]


# Выкидываем рёбра, которые никому не принадлежат, из таблицы:
nobody_edges = df_ref[df_ref["strains"].isnull()].index
print("{n} edges are Nobody".format(n=len(nobody_edges)), '\n')
df_ref = df_ref.loc[set(df_ref.index) - set(nobody_edges)]


# Сплитим список референсов:
df_ref["strains"] = df_ref["strains"].str.split('\t')
df_ref["strains"] = df_ref["strains"].apply(lambda x: [s.rpartition('_')[0] for s in x])
df_ref["strains"] = df_ref["strains"].apply(Counter)


# Считаем копийность каждого ребра:
df_ref["single_copy"] = df_ref["strains"].apply(lambda x: x.most_common(1)[0][1] == 1)


# Считываем длину рёбер:
df_length = pd.read_csv("edges_lengths.tsv", sep='\t', header=None, names=["length"])
df_length = df_length.astype("int")
df_ref = df_ref.join(df_length, how='inner')


# Считываем профили
ref_profile = pd.read_csv("refs/profile.csv", header=None, index_col=0)
for i in range(1, 11):
    ref_profile[i] = ref_profile[i] / ref_profile[i].sum()

desman_profile = pd.read_csv("desman_results/%s_0/desman_freqs.csv" % len(ref_profile),
                             header=None, index_col=0, dtype=float)
desman_profile.index = desman_profile.index.astype(int)


# Ищем соответствие между профилями:
ref_freqs = ref_profile.as_matrix()
ans_error = float("Inf")
ans_permut = None
for cur_permut in permutations(desman_profile.index):
    desman_freqs = desman_profile.loc[cur_permut, :].as_matrix()
    #print(cur_error, cur_permut)
    cur_error = ((ref_freqs - desman_freqs) ** 2).sum()
    if cur_error < ans_error:
        ans_error = cur_error
        ans_permut = cur_permut
print("Error:", ans_error)

def invert_permutation(permutation):
    return [i for i, j in sorted(enumerate(permutation), key=lambda x: x[1])]

strains = list('s' + ref_profile.iloc[invert_permutation(ans_permut), :].index.astype(str))


# Табличка ответов DESMAN
df_desman = pd.read_csv("desman_results/%s_0/gene_assignment_etaS_df.csv" % len(strains), skiprows=1,
                  names=["e_id"] + strains)
df_desman['e_id'] = df_desman['e_id'].str[1:].astype("int")
df_desman = df_desman.set_index('e_id')
df_desman[strains] = df_desman[strains].astype('int')
#df_desman = df_desman.join(df_length, how='inner')
df_desman = df_desman.loc[set(df_ref.index) - set(nobody_edges)]

for cur_s in strains:
    df_ref[cur_s] = df_ref['strains'].apply(lambda x: int(cur_s in x))
    

# Точность DESMAN
df_ref.sort_index(inplace=True)
df_desman.sort_index(inplace=True)

right_answers = (df_ref[strains] == df_desman[strains]).sum(axis=1) == len(strains)
print("Accuracy on all edges: %.2f" % (right_answers.sum() / len(df_ref)))


long_edges = df_ref['length'] > 500
print("Percent of long edges: %.2f" % (long_edges.sum() / len(df_ref)))
print("Accuracy on long edges: %.2f" % ((right_answers & long_edges).sum() / long_edges.sum()))
print()


for cur_s in strains:
    print(cur_s)
    right_answers = df_ref[cur_s] == df_desman[cur_s]
    print("Accuracy on all edges: %.2f" % (right_answers.sum() / len(df_ref)))
    print("Accuracy on long edges: %.2f" % ((right_answers & long_edges).sum() / long_edges.sum()))
    print()


tips_1 = 0
tips_2 = 0
for v in G.nodes:
    if G.in_degree(v) == 0:
        tips_1 += 1
    if G.out_degree(v) == 0:
        tips_2 += 1
print("in_degree == 0, out_degree == 0, all_nodes")
print(tips_1, tips_2, len(G.nodes))
