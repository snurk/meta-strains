
def find_superbubble(G, G_rev, s):
    queue = set([s])
    labels = {}

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
                    return True, t
                else:
                    return False, "cycle"
        if len(queue) == 0:
            break

    return False, "end"
