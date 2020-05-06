import networkx as nx
import pytest
from pseudoflow import hpf


@pytest.fixture()
def G():

    G = nx.DiGraph()
    G.add_edges_from([(0, 1), (1, 2)])

    G[0][1]["const"] = 1
    G[1][2]["const"] = 9.0

    G[0][1]["mult"] = 5
    G[1][2]["mult"] = -3

    return G


#
def test_hpf_nonparametric(G):
    source = 0
    sink = 2
    breakpoints, cuts, info = hpf(G, source, sink, const_cap="const")

    breakpoints == [None]
    cuts == {0: [1], 1: [0], 2: [0]}

    # TODO: Implement test for info! Results are currently incorrect.


def test_hpf_parametric(G):
    from pseudoflow import hpf

    source = 0
    sink = 2
    lambdaRange = [0.0, 2.0]

    breakpoints, cuts, info = hpf(
        G,
        source,
        sink,
        const_cap="const",
        mult_cap="mult",
        lambdaRange=lambdaRange,
        roundNegativeCapacity=False,
    )

    assert breakpoints == [1.0, 2.0]
    assert cuts == {0: [1, 1], 1: [0, 1], 2: [0, 0]}


def test_hpf_with_parametric_sink_arcs():
    from pseudoflow import hpf

    digraph = nx.DiGraph()

    digraph.add_edge("s", 0, const=-20, mult=20)
    digraph.add_edge("s", 1, const=-14, mult=20)
    digraph.add_edge("s", 2, const=-6, mult=20)

    digraph.add_edge(0, "t", const=20, mult=-20)
    digraph.add_edge(1, "t", const=14, mult=-20)
    digraph.add_edge(2, "t", const=6, mult=-20)

    digraph.add_edge(0, 1, const=2, mult=0)
    digraph.add_edge(0, 2, const=1, mult=0)
    digraph.add_edge(2, 1, const=3, mult=0)

    source = "s"
    sink = "t"
    lambdaRange = [0.0, 1.0001]

    breakpoints, cuts, info = hpf(
        digraph,
        source,
        sink,
        const_cap="const",
        mult_cap="mult",
        lambdaRange=lambdaRange,
        roundNegativeCapacity=True,
    )

    assert breakpoints == pytest.approx([0.45, 0.55, 1.0, 1.0001])
    assert cuts == {
        "s": [1, 1, 1, 1],
        0: [0, 0, 0, 1],
        1: [0, 0, 1, 1],
        2: [0, 1, 1, 1],
        "t": [0, 0, 0, 0],
    }


def test_missing_breakpoint():
    G = nx.DiGraph()

    num_deviation = 3.0
    sum_arc_weights = 4.74

    G.add_edge("s", 0, constant=-1.64 / num_deviation, multiplier=2 / num_deviation)
    G.add_edge("s", 3, constant=-0.78 / num_deviation, multiplier=2 / num_deviation)
    G.add_edge("s", 4, constant=-1.02 / num_deviation, multiplier=2 / num_deviation)
    G.add_edge(0, "t", constant=1.64 / num_deviation, multiplier=-2 / num_deviation)
    G.add_edge(3, "t", constant=0.78 / num_deviation, multiplier=-2 / num_deviation)
    G.add_edge(4, "t", constant=1.02 / num_deviation, multiplier=-2 / num_deviation)

    G.add_edge(0, 1, constant=0.88 / sum_arc_weights, multiplier=0)
    G.add_edge(0, 3, constant=0.67 / sum_arc_weights, multiplier=0)
    G.add_edge(1, 0, constant=0.35 / sum_arc_weights, multiplier=0)
    G.add_edge(1, 2, constant=0.24 / sum_arc_weights, multiplier=0)
    G.add_edge(1, 3, constant=0.20 / sum_arc_weights, multiplier=0)
    G.add_edge(1, 4, constant=0.24 / sum_arc_weights, multiplier=0)
    G.add_edge(2, 3, constant=0.92 / sum_arc_weights, multiplier=0)
    G.add_edge(3, 0, constant=0.21 / sum_arc_weights, multiplier=0)
    G.add_edge(3, 1, constant=0.36 / sum_arc_weights, multiplier=0)
    G.add_edge(3, 2, constant=0.12 / sum_arc_weights, multiplier=0)
    G.add_edge(4, 0, constant=0.31 / sum_arc_weights, multiplier=0)
    G.add_edge(4, 3, constant=0.24 / sum_arc_weights, multiplier=0)

    breakpoints, cuts, _ = hpf(
        G,
        "s",
        "t",
        const_cap="constant",
        mult_cap="multiplier",
        lambdaRange=[0.000, 1.0001],
        roundNegativeCapacity=True,
    )
    assert breakpoints == pytest.approx([0.570380, 0.5748101265822786, 1.0001])
    assert cuts == {
        "s": [1, 1, 1],
        "t": [0, 0, 0],
        0: [0, 0, 1],
        1: [0, 0, 1],
        2: [0, 1, 1],
        3: [0, 1, 1],
        4: [0, 0, 1],
    }


def test_problem1():
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

    expected_breakpoints = [0.375049807267, 0.644850603945, 1.0001]
    expected_cuts = {
        "s": [1, 1, 1],
        "t": [0, 0, 0],
        0: [0, 1, 1],
        1: [0, 1, 1],
        2: [0, 1, 1],
        3: [0, 0, 1],
        4: [0, 1, 1],
        5: [0, 1, 1],
        6: [0, 1, 1],
        7: [0, 1, 1],
        8: [0, 0, 1],
        9: [0, 0, 1],
    }

    num_deviation = sum(is_training_dict.values())
    sum_arc_weights = sum(arc_weights.values())

    G = nx.DiGraph()

    for node, is_training in is_training_dict.items():
        if is_training:
            constant = 2 * training_label[node]
            G.add_edge(
                "s",
                node,
                constant=-constant / num_deviation,
                multiplier=2 / num_deviation,
            )
            G.add_edge(
                node,
                "t",
                constant=constant / num_deviation,
                multiplier=-2 / num_deviation,
            )

    for (from_node, to_node), weight in arc_weights.items():
        G.add_edge(from_node, to_node, constant=weight / sum_arc_weights, multiplier=0)

    breakpoints, cuts, _ = hpf(
        G,
        "s",
        "t",
        const_cap="constant",
        mult_cap="multiplier",
        lambdaRange=[0.000, 1.0001],
        roundNegativeCapacity=True,
    )
    assert breakpoints == pytest.approx(expected_breakpoints)
    assert cuts == expected_cuts


def test_problem2():
    is_training_dict = {
        0: True,
        1: True,
        2: True,
        3: False,
        4: True,
        5: False,
        6: False,
        7: False,
        8: False,
        9: False,
    }
    training_label = {
        0: 0.76,
        1: 0.89,
        2: 0.18,
        3: 0.05,
        4: 0.0,
        5: 0.35,
        6: 0.44,
        7: 0.63,
        8: 0.77,
        9: 0.92,
    }

    arc_weights = {
        (0, 1): 0.185,
        (0, 2): 0.38,
        (0, 4): 0.49,
        (0, 5): 0.245,
        (0, 7): 0.25,
        (0, 8): 0.465,
        (1, 0): 0.115,
        (1, 2): 0.41,
        (1, 3): 0.005,
        (1, 5): 0.37,
        (1, 7): 0.08,
        (2, 1): 0.33,
        (2, 3): 0.07,
        (2, 7): 0.23,
        (3, 0): 0.485,
        (3, 1): 0.17,
        (3, 2): 0.23,
        (3, 4): 0.435,
        (3, 6): 0.215,
        (3, 7): 0.18,
        (4, 2): 0.205,
        (4, 5): 0.275,
        (4, 6): 0.295,
        (5, 2): 0.265,
        (5, 4): 0.085,
        (5, 6): 0.015,
        (6, 2): 0.21,
        (6, 3): 0.49,
        (6, 4): 0.05,
        (6, 5): 0.355,
        (6, 7): 0.055,
        (7, 0): 0.185,
        (7, 3): 0.195,
        (7, 4): 0.28,
        (7, 5): 0.195,
        (7, 6): 0.025,
        (7, 8): 0.475,
        (8, 0): 0.31,
        (8, 3): 0.14,
        (9, 2): 0.12,
        (9, 3): 0.16,
        (9, 4): 0.38,
        (9, 6): 0.005,
        (9, 7): 0.26,
    }

    expected_breakpoints = [0.149469, 0.211822, 0.710819, 0.757888, 1.0001]
    expected_cuts = {
        "s": [1, 1, 1, 1, 1],
        "t": [0, 0, 0, 0, 0],
        0: [0, 0, 0, 1, 1],
        1: [0, 0, 0, 0, 1],
        2: [0, 0, 1, 1, 1],
        3: [0, 0, 0, 1, 1],
        4: [0, 1, 1, 1, 1],
        5: [0, 0, 1, 1, 1],
        6: [0, 0, 0, 1, 1],
        7: [0, 0, 0, 1, 1],
        8: [0, 0, 0, 1, 1],
        9: [0, 0, 0, 0, 0],
    }

    num_deviation = sum(is_training_dict.values())
    sum_arc_weights = sum(arc_weights.values())

    G = nx.DiGraph()

    for node, is_training in is_training_dict.items():
        if is_training:
            constant = 2 * training_label[node]
            G.add_edge(
                "s",
                node,
                constant=-constant / num_deviation,
                multiplier=2 / num_deviation,
            )
            G.add_edge(
                node,
                "t",
                constant=constant / num_deviation,
                multiplier=-2 / num_deviation,
            )

    for (from_node, to_node), weight in arc_weights.items():
        G.add_edge(from_node, to_node, constant=weight / sum_arc_weights, multiplier=0)

    breakpoints, cuts, _ = hpf(
        G,
        "s",
        "t",
        const_cap="constant",
        mult_cap="multiplier",
        lambdaRange=[0.000, 1.0001],
        roundNegativeCapacity=True,
    )
    assert breakpoints == pytest.approx(expected_breakpoints, abs=1e-5)
    assert cuts == expected_cuts


def test_problem_seg_fault():
    # only testing for completion.

    is_training_dict = {
        0: False,
        1: False,
        2: True,
        3: False,
        4: False,
        5: False,
        6: True,
        7: True,
        8: False,
        9: True,
    }
    training_label = {
        0: 0.1,
        1: 0.4,
        2: 0.94,
        3: 0.09,
        4: 0.93,
        5: 0.47,
        6: 0.45,
        7: 0.44,
        8: 0.27,
        9: 0.14,
    }

    arc_weights = {
        (0, 2): 0.345,
        (0, 3): 0.065,
        (0, 4): 0.155,
        (0, 7): 0.13,
        (0, 8): 0.08,
        (1, 0): 0.335,
        (1, 2): 0.255,
        (1, 3): 0.455,
        (1, 5): 0.065,
        (1, 6): 0.22,
        (1, 7): 0.255,
        (1, 9): 0.495,
        (2, 0): 0.18,
        (2, 4): 0.035,
        (2, 6): 0.36,
        (2, 8): 0.14,
        (2, 9): 0.24,
        (3, 1): 0.465,
        (3, 2): 0.49,
        (3, 4): 0.28,
        (3, 6): 0.265,
        (3, 7): 0.19,
        (3, 9): 0.37,
        (4, 0): 0.485,
        (4, 2): 0.045,
        (4, 7): 0.27,
        (4, 8): 0.29,
        (5, 0): 0.22,
        (5, 4): 0.17,
        (5, 7): 0.15,
        (6, 0): 0.23,
        (6, 2): 0.195,
        (6, 3): 0.175,
        (6, 4): 0.44,
        (6, 5): 0.48,
        (6, 7): 0.2,
        (6, 8): 0.065,
        (7, 0): 0.19,
        (7, 4): 0.35,
        (7, 5): 0.39,
        (7, 8): 0.225,
        (8, 1): 0.39,
        (8, 3): 0.34,
        (8, 4): 0.33,
        (8, 5): 0.21,
        (8, 6): 0.215,
        (8, 9): 0.035,
        (9, 1): 0.16,
        (9, 4): 0.07,
        (9, 5): 0.065,
        (9, 7): 0.11,
    }

    num_deviation = sum(is_training_dict.values())
    sum_arc_weights = sum(arc_weights.values())

    G = nx.DiGraph()

    for node, is_training in is_training_dict.items():
        if is_training:
            constant = 2 * training_label[node]
            G.add_edge(
                "s",
                node,
                constant=-constant / num_deviation,
                multiplier=2 / num_deviation,
            )
            G.add_edge(
                node,
                "t",
                constant=constant / num_deviation,
                multiplier=-2 / num_deviation,
            )

    for (from_node, to_node), weight in arc_weights.items():
        G.add_edge(from_node, to_node, constant=weight / sum_arc_weights, multiplier=0)

    breakpoints, cuts, _ = hpf(
        G,
        "s",
        "t",
        const_cap="constant",
        mult_cap="multiplier",
        lambdaRange=[0.000, 1.0001],
        roundNegativeCapacity=True,
    )


def test_sink_adjacent_arcs_with_positive_multiplier_raises_valueerror():
    G = nx.DiGraph()
    G.add_nodes_from(["source", "sink", 0, 1, 2])

    G.add_edge(0, 1, weight=float("inf"), multiplier=0)
    G.add_edge(0, 2, weight=float("inf"), multiplier=0)
    G.add_edge(1, 2, weight=float("inf"), multiplier=0)

    G.add_edge(0, "sink", weight=3, multiplier=0)
    G.add_edge("source", 1, weight=6, multiplier=0)
    G.add_edge(2, "sink", weight=4, multiplier=2)

    with pytest.raises(ValueError):
        hpf(
            G,
            "source",
            "sink",
            const_cap="weight",
            mult_cap="multiplier",
            lambdaRange=[0, 5],
            roundNegativeCapacity=False,
        )


def test_source_adjacent_arcs_with_negative_multiplier_raises_valueerror():
    G = nx.DiGraph()
    G.add_nodes_from(["source", "sink", 0, 1, 2])

    G.add_edge(0, 1, weight=float("inf"), multiplier=0)
    G.add_edge(0, 2, weight=float("inf"), multiplier=0)
    G.add_edge(1, 2, weight=float("inf"), multiplier=0)

    G.add_edge(0, "sink", weight=3, multiplier=0)
    G.add_edge("source", 1, weight=6, multiplier=-2)
    G.add_edge(2, "sink", weight=4, multiplier=0)

    with pytest.raises(ValueError):
        hpf(
            G,
            "source",
            "sink",
            const_cap="weight",
            mult_cap="multiplier",
            lambdaRange=[0, 5],
            roundNegativeCapacity=False,
        )
