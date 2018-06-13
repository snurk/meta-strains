def print_accuracy(cur_s, cur_ans):
    cur_true = cur_ans[cur_s] == 1
    G_sub = G.subgraph(to_double_format(cur_ans[cur_true].index))
    print("components in DESMAN  subgraph:", nx.number_weakly_connected_components(G_sub))

    right_answers = df_ref[cur_s] == cur_ans[cur_s]
    print("Accuracy on all edges: %.2f" % (right_answers.sum() / len(df_ref)))


def print_stats(cur_s, cur_df, ans):
    real_true = df_ref[cur_s] == 1
    desman_true = cur_df[cur_s] == 1

    print("\nold FN")
    print_nodes(df_ref[real_true & ~desman_true].index)

    print("\nold FP")
    print_nodes(df_ref[~real_true & desman_true].index)

    print("\nremained FN:")
    print_nodes(set(df_ref[real_true & ~desman_true].index) - ans)

    # paths = list(nx.all_simple_paths(G, source=s, target=t))
    # bubble_sizes.append(len(paths))
    # print(t)
    # bubbles_s_t.append((s, t, rev_flag))
    # print("!", len(list(paths)), s, t)
    # print('!', set(x for lst in list(paths) for x in lst[1:-1]))