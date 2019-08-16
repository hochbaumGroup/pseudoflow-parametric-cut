import pytest
import networkx as nx
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

    expected_breakpoints = [0.375049807267, 0.532261662979, 0.644850603945, 1.0001]
    expected_cuts = {
        "s": [1, 1, 1, 1],
        "t": [0, 0, 0, 0],
        0: [0, 1, 1, 1],
        1: [0, 1, 1, 1],
        2: [0, 1, 1, 1],
        3: [0, 0, 0, 1],
        4: [0, 1, 1, 1],
        5: [0, 1, 1, 1],
        6: [0, 1, 1, 1],
        7: [0, 1, 1, 1],
        8: [0, 0, 0, 1],
        9: [0, 0, 1, 1],
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
