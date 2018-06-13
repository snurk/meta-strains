# coding: utf-8
import itertools
import numpy as np

from read_files import *
from superbubbles import find_superbubble

import multiprocessing
import sys

import functools

def print_stats(df_ref, df_ans, G, title=None, dataset=""):
    print()
    if title:
        print(title)

    precision_sum = 0
    recall_sum = 0
    gf_sum = 0
    trash_len_sum = 0

    strains = df_ans.columns
    k = 126
    if dataset == "infant_gut":
        strains = ["s1", "s3"]
        k = 76

    for cur_s in strains:
        real_true = df_ref[cur_s] == 1
        predicted_true = df_ans[cur_s] == 1

        TP_df = df_ref[real_true & predicted_true]
        TP = len(TP_df)
        #print_nodes(TP_df.index)

        FP_df  = df_ref[~real_true & predicted_true]
        FP = len(FP_df)

        TN_df  = df_ref[~real_true & ~predicted_true]
        TN = len(TN_df)

        FN_df  = df_ref[real_true & ~predicted_true]
        FN = len(FN_df)

        precision_sum += TP / (TP+FP+1e-10)
        recall_sum += TP / (TP+FN+1e-10)

        ref_len = df_ref.loc[real_true, "length"].sum() - k * (real_true.sum() - 1)
        TP_len = TP_df["length"].sum() - k * (TP - 1)
        genome_fraction = TP_len / ref_len
        gf_sum += genome_fraction
        
        FP_len = FP_df["length"].sum() - k * (FP - 1)
        trash_len_sum += FP_len

        if dataset == "infant_gut":
            print("________", cur_s)
            print("TP={}  TN={}  FP={}  FN={}".format(TP, TN, FP, FN))
            print("precision: {0:.2f}".format(TP / (TP+FP+1e-10)))
            print("   recall: {0:.2f}".format(TP / (TP+FN+1e-10)))
            print("       GF: {0:.2f}".format(genome_fraction))
            print("FP  len: {}".format(FP_len))


        #print("ref len: {}".format(ref_len))
        #print("TP  len: {}".format(TP_len))
        #print("FP  len: {}".format(FP_len))

        #G_sub = G.subgraph(to_double_format(df_ans[predicted_true].index))
        #print("components in DESMAN  subgraph:", nx.number_weakly_connected_components(G_sub))

    if dataset == "infant_gut":
        return

    precision_mean = precision_sum / len(strains)
    recall_mean = recall_sum / len(strains)
    gf_mean = gf_sum / len(strains)
    trash_len_mean = trash_len_sum / len(strains)
    print("precision: {0:.2f}".format(precision_mean))
    print("   recall: {0:.2f}".format(recall_mean))
    print("       GF: {0:.2f}".format(gf_mean))
    print("FP len: {0:.0f}".format(trash_len_mean))


def correct_cutpoints(G, df_ans):
    G_rev = G.reverse()

    df_corrected_ans = df_ans.copy()
    bubbles_s_t = []

    for cur_s in df_ans.columns:
        selected_nodes = to_double_format(df_ans[df_ans[cur_s] == 1].index)
        selected_nodes = [s for s in selected_nodes if s[-1] == '+']

        visited = dict.fromkeys(selected_nodes, False)
        ans = []

        for node in selected_nodes:
            if visited[node] or visited.get(rev(node), False):
                continue

            visited[node] = True

            for graph, rev_flag in zip([G, G_rev], [False, True]):

                if not rev_flag:
                    rev_graph = G_rev
                else:
                    rev_graph = G

                cur_node = node
                has_bubble, t_bubb_end = find_superbubble(graph, rev_graph, node)
                unique_continue, t_unique_cont = graph.out_degree(cur_node) == 1, list(graph.neighbors(cur_node))
                while has_bubble or unique_continue:
                    if has_bubble:
                        bubbles_s_t.append((cur_node, t_bubb_end, rev_flag))
                        t = t_bubb_end
                    else:  # unique_continue
                        t = t_unique_cont[0]

                    if visited.get(t, False):
                        break
                    visited[t] = True
                    ans.append(t)
                    cur_node = t

                    has_bubble, t_bubb_end = find_superbubble(graph, rev_graph, cur_node)
                    unique_continue, t_unique_cont = graph.out_degree(cur_node) == 1, list(graph.neighbors(cur_node))

        ans = to_single_format(ans)
        df_corrected_ans.loc[ans, cur_s] = 1

    return df_corrected_ans, bubbles_s_t


def correct_bubbles(G, df_ans, bubbles_s_t, desman_profile):
    n_samples = desman_profile.shape[1]
    strains = df_ans.columns

    assignments = {}

    #print('\n', "number of bubbles:", len(bubbles_s_t), '\n')
    iter_step = 0
    for s, t, rev_flag in bubbles_s_t:

        iter_step += 1
        if iter_step % 200 == 0:
            print(iter_step)
        #if iter_step == 1000:
        #    break

        if rev_flag:
            s, t = t, s
        paths = list(nx.all_simple_paths(G, source=s, target=t))
        paths = [p[1:-1] for p in paths]
        inner_nodes = set(x for lst in list(paths) for x in lst)
        cur_strains_bool = df_ans.loc[to_my_int(s)]
        cur_strains = [s for s in strains if cur_strains_bool[s]]

        min_error = float("Inf")
        min_variant = None

        cur_profile = desman_profile.loc[cur_strains]
        cur_profile /= cur_profile.sum()
        paths_s = [paths for i in range(len(cur_strains))]

        if len(paths) ** len(cur_strains) > 3000:
            print('!', len(paths), len(paths) ** len(cur_strains))
            continue

        mean_cov = G.node[s]["cov"] * G.node[s]["length"] + G.node[t]["cov"] * G.node[t]["length"] + 1e-10
        mean_cov /= G.node[s]["length"] + G.node[t]["length"] #- 2 * 126

        for variant in itertools.product(*paths_s):

            nodes_dict = {}
            for cur_node in inner_nodes:
                nodes_dict[cur_node] = {}
                nodes_dict[cur_node]["strains"] = []
                nodes_dict[cur_node]["real_cov"] = G.node[cur_node]["cov"]
                nodes_dict[cur_node]["estimated_prof"] = np.zeros(n_samples)
                nodes_dict[cur_node]["real_prof"] = np.zeros(n_samples)

            for i in range(len(cur_strains)):
                s_i = cur_strains[i]
                for cur_node in variant[i]:
                    nodes_dict[cur_node]["strains"].append(s_i)
                    nodes_dict[cur_node]["estimated_prof"] += cur_profile.loc[s_i]

            cur_error = np.zeros(n_samples)
            for node in nodes_dict:
                nodes_dict[node]["real_prof"] = nodes_dict[node]["real_cov"] / mean_cov
                node_error = np.linalg.norm(nodes_dict[node]["estimated_prof"] - nodes_dict[node]["real_prof"]) ** 2
                cur_error += node_error * np.log(G.node[node]["length"])
                if cur_error.sum() >= min_error:
                    break

            cur_error = cur_error.sum()

            if cur_error < min_error:
                min_error = cur_error
                for cur_node in inner_nodes:
                    nodes_dict[cur_node] = nodes_dict[cur_node]["strains"]
                min_variant = nodes_dict

        assignments.update(min_variant)

    return assignments


def merge_bubbles_answers(df_ans, assignments):
    df_corrected_ans = df_ans.copy()
    for node, cur_strains in assignments.items():
        df_corrected_ans.loc[to_my_int(node), cur_strains] = 1
    return df_corrected_ans


def main():
    #dataset = "g3_r1"
    dataset = sys.argv[1]
    print(dataset)

    G = read_graph(dataset)

    df_ref, df_desman, desman_profile = read_answers(G, dataset)
    print_stats(df_ref, df_desman, G, "*****Initial", dataset)

    df_corrected_cutpoints, bubbles_s_t = correct_cutpoints(G, df_desman)
    print_stats(df_ref, df_corrected_cutpoints, G, "*****Continues", dataset)

    chunks = [bubbles_s_t[x:x + 100] for x in range(0, len(bubbles_s_t), 100)]
    bubbles_assignments = {}

    pool = multiprocessing.Pool(int(sys.argv[2]))
    n = len(chunks)
    args = zip([G] * n, [df_corrected_cutpoints] * n, chunks, [desman_profile] * n)
    assignments = pool.starmap(correct_bubbles, args)
    #assignments = pool.map(functools.partial(correct_bubbles, G=G, df_corrected_cutpoints=df_corrected_cutpoints, desman_profile=desman_profile), chunks)
    for ans in assignments:
       bubbles_assignments.update(ans)

    df_corrected_bubbles = merge_bubbles_answers(df_corrected_cutpoints, bubbles_assignments)
    print_stats(df_ref, df_corrected_bubbles, G, "*****Bubbles", dataset)


if __name__ == '__main__':
    main()
