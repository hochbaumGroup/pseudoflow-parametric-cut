import pytest
import networkx as nx


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
    from pseudoflow import hpf

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

    print(breakpoints)
    print(cuts)
    assert breakpoints == pytest.approx([0.45, 0.55, 1.0, 1.0001])
    assert cuts == {
        "s": [1, 1, 1, 1],
        0: [0, 0, 0, 1],
        1: [0, 0, 1, 1],
        2: [0, 1, 1, 1],
        "t": [0, 0, 0, 0],
    }
