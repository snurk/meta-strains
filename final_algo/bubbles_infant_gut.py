# coding: utf-8


from read_files import *


def find_superbubbles_ends(G, G_rev, s):
    queue = set([s])
    labels = {}
    n = 0

    while True:

        v = queue.pop()
        labels[v] = "visited"

        v_children = list(G.neighbors(v))

        if len(v_children) == 0:
            return False, "tip"

        for u in v_children:
            if u == s:
                return False, "cycle"

            labels[u] = "seen"

            u_parents = list(G_rev.neighbors(u))

            if all(labels.get(parent, "") == "visited" for parent in u_parents):
                queue.add(u)

        if len(queue) == 1:
            t = list(queue)[0]
            if "seen" not in [labels[k] for k in set(labels.keys()) - set([t])]:

                if not G.has_edge(s, t):
                    paths = nx.all_simple_paths(G, source=s, target=t)
                    print("!", len(list(paths)), s, t)
                    return True, t
                else:
                    return False, "cycle"
        if len(queue) == 0:
            break

    return False, "end"


# In[ ]:

G = read_graph("infant_gut_graph.gfa")

#df_ref, df_desman, strains = read_answers(G)


G_rev = G.reverse()

bubble_ends = set()
for node in G.nodes:
    has_bubble, bubble_end = find_superbubbles_ends(G, G_rev, node)
    if has_bubble:
        bubble_ends.add(node[:-1])
        bubble_ends.add(bubble_end[:-1])
#print_nodes(bubble_ends)


simple_repeats = set()
for node in G.nodes:
    if G.in_degree(node) == 1 and G.out_degree(node) == 1:
        if list(G.predecessors(node))[0] == list(G.successors(node))[0]:
            simple_repeats.add(node[:-1])
print_nodes(simple_repeats)



#df_desman_corrected = df_desman.copy()

# for cur_s in strains:
#
#     print('\n\n', "_________________________", cur_s, '\n')
#
#     ref_true = df_ref[cur_s] == 1
#     G_sub = G.subgraph(to_double_format(df_ref[ref_true].index))
#     print("components in reference subgraph:", nx.number_weakly_connected_components(G_sub))
#
#     desman_true = df_desman[cur_s] == 1
#     G_sub = G.subgraph(to_double_format(df_desman[desman_true].index))
#     print("components in   DESMAN  subgraph:", nx.number_weakly_connected_components(G_sub))
#
#     long_1500 = df_ref['length'] >= 1500
#
#     selected_nodes = to_double_format(df_desman[long_1500 & desman_true].index)
#     selected_nodes = [s for s in selected_nodes if (s[-1] == '+' and len(s) < 12)]
#
#     visited = dict.fromkeys(selected_nodes, False)
#     bubbles = []
#     ans = []
#
#     print("selected nodes", len(selected_nodes))
#
#     for node in selected_nodes:
#         if visited[node] or visited.get(rev(node), False):
#             continue
#
#         visited[node] = True
#
#         for graph in [G, G_rev]:
#
#             if graph is G:
#                 rev_graph = G_rev
#             else:
#                 rev_graph = G
#
#             has_bubble, bubble_end = find_superbubble(graph, rev_graph, node)
#             while has_bubble:
#                 if visited.get(bubble_end, False):
#                     break
#                 visited[bubble_end] = True
#                 if df_desman.loc[to_my_int(bubble_end), cur_s] == 0:
#                     ans.append(bubble_end)
#                 has_bubble, bubble_end = find_superbubble(graph, rev_graph, bubble_end)
#
#     ans = to_single_format(ans)
#     df_desman_corrected.loc[ans, cur_s] = 1
#
#     right_answers_old = df_ref[cur_s] == df_desman[cur_s]
#     print("Old accuracy on all edges: %.2f" % (right_answers_old.sum() / len(df_ref)))
#
#     right_answers_new = df_ref[cur_s] == df_desman_corrected[cur_s]
#     print("New accuracy on all edges: %.2f" % (right_answers_new.sum() / len(df_ref)))
#
#     G_sub = G.subgraph(to_double_format(df_desman_corrected[df_desman_corrected[cur_s] == 1].index))
#     print("components in new DESMAN subgraph:", nx.number_weakly_connected_components(G_sub))
#
#     real_true = df_ref[cur_s] == 1
#     desman_true = df_desman[cur_s] == 1
#
#     print("\nold FN")
#     print_nodes(df_ref[real_true & ~desman_true].index)
#
#     print("\nold FP")
#     print_nodes(df_ref[~real_true & desman_true].index)
#
#     print("\nremained FN:")
#     print_nodes(set(df_ref[real_true & ~desman_true].index) - ans)
#
#     print("\nnew FP:")
#     print_nodes(ans - set(df_ref[real_true].index))

# In[ ]:
