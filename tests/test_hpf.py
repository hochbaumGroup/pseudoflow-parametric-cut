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

    G.add_edge("s", 0, constant=-1.64 / 3.0, multiplier=2 / 3.0)
    G.add_edge("s", 3, constant=-0.78 / 3.0, multiplier=2 / 3.0)
    G.add_edge("s", 4, constant=-1.02 / 3.0, multiplier=2 / 3.0)
    G.add_edge(0, "t", constant=1.64 / 3.0, multiplier=-2 / 3.0)
    G.add_edge(3, "t", constant=0.78 / 3.0, multiplier=-2 / 3.0)
    G.add_edge(4, "t", constant=1.02 / 3.0, multiplier=-2 / 3.0)

    G.add_edge(0, 1, constant=0.88 / 4.74, multiplier=0)
    G.add_edge(0, 3, constant=0.67 / 4.74, multiplier=0)
    G.add_edge(1, 0, constant=0.35 / 4.74, multiplier=0)
    G.add_edge(1, 2, constant=0.24 / 4.74, multiplier=0)
    G.add_edge(1, 3, constant=0.20 / 4.74, multiplier=0)
    G.add_edge(1, 4, constant=0.24 / 4.74, multiplier=0)
    G.add_edge(2, 3, constant=0.92 / 4.74, multiplier=0)
    G.add_edge(3, 0, constant=0.21 / 4.74, multiplier=0)
    G.add_edge(3, 1, constant=0.36 / 4.74, multiplier=0)
    G.add_edge(3, 2, constant=0.12 / 4.74, multiplier=0)
    G.add_edge(4, 0, constant=0.31 / 4.74, multiplier=0)
    G.add_edge(4, 3, constant=0.24 / 4.74, multiplier=0)

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
