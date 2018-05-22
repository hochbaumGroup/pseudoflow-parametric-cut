import pytest


@pytest.fixture()
def G():
    from networkx import DiGraph

    G = DiGraph()
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
    breakpoints, cuts, info = hpf(
        G, source, sink, const_cap="const")

    breakpoints == [None, ]
    cuts == {0: [1, ], 1: [0, ], 2: [0, ]}

    # TODO: Implement test for info! Results are currently incorrect.


def test_hpf_parametric(G):
    from pseudoflow import hpf

    source = 0
    sink = 2
    lambdaRange = [0., 2.]

    breakpoints, cuts, info = hpf(
        G, source, sink, const_cap="const", mult_cap="mult",
        lambdaRange=lambdaRange, roundNegativeCapacity=False)

    assert breakpoints == [1., 2.]
    assert cuts == {0: [1, 1], 1: [0, 1], 2: [0, 0]}
