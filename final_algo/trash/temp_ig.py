import pandas as pd
import csv
import networkx as nx
from graph_functions import *

from read_files import read_graph


def read_answers_1(G, dataset_name="example"):
    df_ref = pd.read_csv("data/{}/refs_edges_0.txt".format(dataset_name), header=None, names=["e"])

    df_ref = df_ref["e"].str.split('\t', 1, expand=True)
    df_ref.columns = ["e_id", "strains"]
    df_ref = df_ref.set_index("e_id")

    # remain only nodes from largest component
    df_ref = df_ref.loc[set([to_my_int(node) for node in G.nodes])]


    # DESMAN gene_assignment answer
    df_desman = pd.read_csv("data/{}/gene_assignment_etaS_df.csv".format(dataset_name), index_col=0)
    df_desman.index = df_desman.index.str[1:]

    df_desman = df_desman.reindex(df_ref.index)

    #df_desman.loc[set(df_ref.index) - set(df_desman.index)] = 0

    df_desman = df_desman.fillna(0)

    df_desman.index = 'e' + df_desman.index
    df_desman.to_csv("data/{}/gene_assignment_etaS_df_0.csv".format(dataset_name))



G = read_graph("infant_gut")
read_answers_1(G, "infant_gut")
