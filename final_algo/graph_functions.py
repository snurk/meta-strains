def print_nodes(nodes_list):
    print(','.join([str(node) for node in nodes_list]))


def to_single_format(nodes_list):
    return set([to_my_int(node) for node in nodes_list])


def to_double_format(nodes_list):
    nodes_1 = [str(node) + "+" for node in nodes_list]
    nodes_2 = [str(node) + "-" for node in nodes_list]
    return nodes_1 + nodes_2


def to_my_int(node_name):
    node_name = node_name[:-1]
    node_name = "".join(node_name.split('_'))
    return node_name


def rev(node_name):
    direction = node_name[-1]
    if direction == "-":
        rev_direction = "+"
    else:
        rev_direction = "-"
    return node_name[:-1] + rev_direction
