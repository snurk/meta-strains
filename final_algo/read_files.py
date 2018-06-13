import pandas as pd
from collections import Counter
from itertools import permutations
import networkx as nx
from graph_functions import *


def read_graph(dataset_name="example"):
    G = nx.DiGraph()

    df_cov = pd.read_csv("data/{}/edge_profiles_0.txt".format(dataset_name),
                         sep=' ', index_col=0, header=None)
    df_cov.index = df_cov.index.astype(str)

    with open("data/{}/merged_graph_0.gfa".format(dataset_name)) as f:
        for line in f:
            line = line.split()
            line_type = line[0]

            # S 238024 ACCAATTAT KC:i:37210
            if line_type == "S":
                v_name = "".join(line[1].split('_'))
                v_length = len(line[2])
                G.add_node(v_name + "+", length=v_length, cov=df_cov.loc[v_name].values)
                G.add_node(v_name + "-", length=v_length, cov=df_cov.loc[v_name].values)

            # L 238322 + 19590 - 55M
            if line_type == "L":
                v1 = "".join(line[1].split('_')) + line[2]
                v2 = "".join(line[3].split('_')) + line[4]
                G.add_edge(v1, v2)
                G.add_edge(rev(v2), rev(v1))

    # remain only largest component
    Gc = max(nx.weakly_connected_component_subgraphs(G), key=len)
    return Gc


def invert_permutation(permutation):
    return [i for i, j in sorted(enumerate(permutation), key=lambda x: x[1])]


def read_answers(G, dataset_name="example"):
    df_ref = pd.read_csv("data/{}/refs_edges_0.txt".format(dataset_name), header=None, names=["e"])

    df_ref = df_ref["e"].str.split('\t', 1, expand=True)
    df_ref.columns = ["e_id", "strains"]
    df_ref = df_ref.set_index("e_id")

    # remain only nodes from largest component
    df_ref = df_ref.loc[set([to_my_int(node) for node in G.nodes])]

    # search unaligned nodes
    nobody_nodes = df_ref[df_ref["strains"].isnull()].index
    print("{n} nodes are Nobody".format(n=len(nobody_nodes)), '\n')

    df_ref.loc[nobody_nodes, "strains"] = "Nobody_0"

    # parsing references alignment
    df_ref["strains"] = df_ref["strains"].str.split('\t')
    df_ref["strains"] = df_ref["strains"].apply(lambda x: [s.rpartition('_')[0] for s in x])
    df_ref["strains"] = df_ref["strains"].apply(Counter)

    df_ref["single_copy"] = df_ref["strains"].apply(lambda x: x.most_common(1)[0][1] == 1)

    # node lengths calculation
    dict_length = {}

    for k, v in nx.get_node_attributes(G, 'length').items():
        if k[-1] == '+':
            dict_length[to_my_int(k)] = v

    df_length = pd.Series(dict_length, name='length')
    df_length.index.name = 'e_id'
    df_ref = df_ref.join(df_length, how='inner')

    # compare DESMAN and real profiles
    ref_profile = pd.read_csv("data/{}/profile.csv".format(dataset_name), header=None, index_col=0)
    for i in range(1, 11):
        ref_profile[i] = ref_profile[i] / ref_profile[i].sum()
    desman_profile = pd.read_csv("data/{}/desman_freqs.csv".format(dataset_name), header=None, index_col=0, dtype=float)
    desman_profile.index = desman_profile.index.astype(int)
    ans_error, ans_permut = float("Inf"), None
    for cur_permut in permutations(desman_profile.index):
        desman_freqs = desman_profile.loc[cur_permut, :].values
        cur_error = ((ref_profile.values - desman_freqs) ** 2).sum()
        if cur_error < ans_error:
            ans_error, ans_permut = cur_error, cur_permut
    print("Error:%.3f" % ans_error)

    strains = list('s' + ref_profile.iloc[invert_permutation(ans_permut), :].index.astype(str))
    desman_profile.index = strains

    # DESMAN gene_assignment answer
    df_desman = pd.read_csv("data/{}/gene_assignment_etaS_df_0.csv".format(dataset_name), skiprows=1,
                            names=["e_id"] + strains)
    df_desman['e_id'] = df_desman['e_id'].str[1:]
    df_desman = df_desman.set_index('e_id')
    df_desman[strains] = df_desman[strains].astype('int')

    df_desman = df_desman.loc[df_ref.index]

    df_desman[df_ref['length'] < 1000] = 0

    # references alignment answer
    for cur_s in strains:
        df_ref[cur_s] = df_ref['strains'].apply(lambda x: int(cur_s in x))

    df_ref.sort_index(inplace=True)
    df_desman.sort_index(inplace=True)

    # initial accuracy of DESMAN
    right_answers = (df_ref[strains] == df_desman[strains]).sum(axis=1) == len(strains)
    print("Accuracy on all edges: %.2f" % (right_answers.sum() / len(df_ref)))
    long_nodes = df_ref['length'] >= 1000
    print("Percent of long nodes: %.2f" % (long_nodes.sum() / len(df_ref)))
    print("Accuracy on long nodes: %.2f" % ((right_answers & long_nodes).sum() / long_nodes.sum()))

    return df_ref, df_desman, desman_profile
