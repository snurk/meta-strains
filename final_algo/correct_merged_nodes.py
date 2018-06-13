import pandas as pd
import csv
import networkx as nx
from graph_functions import *


def read_graph_0(dataset_name="example"):
    G = nx.DiGraph()

    with open("data/{}/merged_graph.gfa".format(dataset_name)) as f:
        for line in f:
            line = line.split()
            line_type = line[0]

            # S 238024 ACCAATTAT KC:i:37210
            if line_type == "S":
                v_name = line[1]
                v_length = len(line[2])
                G.add_node(v_name + "+", length=v_length)
                G.add_node(v_name + "-", length=v_length)

            # L 238322 + 19590 - 55M
            if line_type == "L":
                v1 = line[1] + line[2]
                v2 = line[3] + line[4]
                G.add_edge(v1, v2)
                G.add_edge(rev(v2), rev(v1))

    # remain only largest component
    new_G = nx.DiGraph()
    for g in nx.weakly_connected_component_subgraphs(G):
        if new_G.number_of_nodes() < g.number_of_nodes():
            new_G = g.copy()
    return new_G


def read_answers_0(G, dataset_name="example"):
    df_ref = pd.read_csv("data/{}/refs_edges.txt".format(dataset_name), header=None, names=["e"])

    df_ref = df_ref["e"].str.split('\t', 1, expand=True)
    df_ref.columns = ["e_id", "strains"]
    df_ref = df_ref.set_index("e_id")

    # assign merged nodes
    for node in G.nodes:
        if '_' in node and '-' in node:
            parts = node[:-1].split('_')
            cur_strains = [df_ref.loc[p, "strains"] for p in parts]
            if not len(set(cur_strains)) == 1:
                print(cur_strains)
            cur_strains = ['' if v is None else v for v in cur_strains]

            df_ref.loc[to_my_int(node), "strains"] = max(cur_strains, key=len)
    df_ref.to_csv("data/{}/refs_edges_0.txt".format(dataset_name),
                  sep='\t', header=False, escapechar="\t", quoting=csv.QUOTE_NONE)

    # DESMAN gene_assignment answer
    df_desman = pd.read_csv("data/{}/gene_assignment_etaS_df.csv".format(dataset_name), index_col=0)
    df_desman.index = df_desman.index.str[1:]

    # assign merged edges
    for node in G.nodes:
        if '_' in node and '-' in node:
            parts = node[:-1].split('_')
            cur_strains = df_desman.loc[[p for p in parts], :]
            df_desman.loc[to_my_int(node)] = cur_strains.min()

            if not (cur_strains.min() == cur_strains.max()).all():
                print(cur_strains)

    df_desman.index = 'e' + df_desman.index
    df_desman.to_csv("data/{}/gene_assignment_etaS_df_0.csv".format(dataset_name))

    # coverage
    df_cov = pd.read_csv("data/{}/edge_profiles.txt".format(dataset_name),
                         sep=' ', index_col=0)
    df_cov.index = df_cov.index.astype(str)
    df_cov.drop(df_cov.columns[len(df_cov.columns) - 1], axis=1, inplace=True)
    # assign merged nodes
    for node in G.nodes:
        if '_' in node and '-' in node:
            parts = node[:-1].split('_')
            cur_strains = df_cov.loc[[p for p in parts], :]
            df_cov.loc[to_my_int(node)] = cur_strains.mean()
    df_cov.to_csv("data/{}/edge_profiles_0.txt".format(dataset_name), sep=' ')


for dataset in ['g2_r1',  'g2_r2',  'g3_r1',  'g3_r2',  'g4_r1',  'g4_r2',  'g5_r1',  'g5_r2']:
    #dataset = "g5_r2"
    G = read_graph_0(dataset)
    read_answers_0(G, dataset)
