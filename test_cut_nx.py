import networkx as nx

is_training_dict = {
    0: True,
    1: True,
    2: True,
    3: True,
    4: True,
    5: True,
    6: False,
    7: True,
    8: True,
    9: False,
}
training_label = {
    0: 0.22,
    1: 0.45,
    2: 0.16,
    3: 0.97,
    4: 0.31,
    5: 0.14,
    6: 0.53,
    7: 0.36,
    8: 0.93,
    9: 0.47,
}

arc_weights = {
    (0, 2): 0.355,
    (0, 7): 0.005,
    (1, 0): 0.315,
    (1, 2): 0.06,
    (1, 4): 0.07,
    (1, 5): 0.005,
    (1, 6): 0.145,
    (1, 7): 0.365,
    (2, 0): 0.455,
    (2, 1): 0.19,
    (2, 3): 0.045,
    (2, 4): 0.49,
    (2, 6): 0.44,
    (2, 7): 0.14,
    (2, 8): 0.23,
    (2, 9): 0.11,
    (3, 1): 0.34,
    (3, 2): 0.195,
    (3, 4): 0.075,
    (3, 5): 0.065,
    (3, 8): 0.395,
    (3, 9): 0.255,
    (4, 0): 0.03,
    (4, 1): 0.04,
    (4, 2): 0.295,
    (4, 5): 0.015,
    (5, 2): 0.04,
    (5, 3): 0.355,
    (5, 6): 0.05,
    (5, 9): 0.05,
    (6, 3): 0.325,
    (6, 8): 0.085,
    (7, 0): 0.295,
    (8, 2): 0.315,
    (8, 3): 0.23,
    (8, 4): 0.35,
    (8, 5): 0.005,
    (8, 6): 0.12,
    (8, 9): 0.25,
    (9, 1): 0.11,
    (9, 8): 0.16,
}


num_deviation = sum(is_training_dict.values())
sum_arc_weights = sum(arc_weights.values())
lambda_value = 0.4426  # 0.4424

G = nx.DiGraph()

for node, is_training in is_training_dict.items():
    if is_training:
        weight = 2 * (lambda_value - training_label[node]) / num_deviation
        G.add_edge("s", node, weight=max(weight, 0))
        G.add_edge(node, "t", weight=max(-weight, 0))

for (from_node, to_node), weight in arc_weights.items():
    G.add_edge(from_node, to_node, weight=weight / sum_arc_weights)

cut_value, partition = nx.minimum_cut(G, "s", "t", capacity="weight")
reachable, non_reachable = partition
print(cut_value, reachable, non_reachable)
