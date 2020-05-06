from ctypes import c_int, c_double, cast, byref, POINTER, cdll
from six.moves import xrange
import os

PATH = os.path.dirname(__file__)
libhpf = cdll.LoadLibrary(os.path.join(PATH, os.pardir, "libhpf.so"))


def _c_arr(c_type, size, init):
    x = c_type * size
    return x(*init)


def _get_arcmatrix(G, const_cap, mult_cap, source, sink):
    nNodes = G.number_of_nodes()
    nArcs = G.number_of_edges()

    nodeNames = []
    nodeDict = {}
    linearArcMatrix = []

    nodeDict[source] = 0
    nodeNames = [source]

    for node in G.nodes():
        if node not in {source, sink}:
            nodeDict[node] = len(nodeNames)
            nodeNames.append(node)

    nodeDict[sink] = len(nodeNames)
    nodeNames.append(sink)

    for fromNode, toNode, data in G.edges(data=True):
        linearArcMatrix += [
            nodeDict[fromNode],
            nodeDict[toNode],
            data[const_cap],
            data[mult_cap] if mult_cap else 0,
        ]

    return (nodeNames, nodeDict, map(lambda x: float(x), linearArcMatrix))


def _create_c_input(
    G, nodeDict, source, sink, arcMatrix, lambdaRange, roundNegativeCapacity
):
    nNodes = G.number_of_nodes()
    nArcs = G.number_of_edges()
    c_numNodes = c_int(nNodes)
    c_numArcs = c_int(nArcs)
    c_source = c_int(nodeDict[source])
    c_sink = c_int(nodeDict[sink])
    c_arcMatrix = _c_arr(c_double, nArcs * 4, arcMatrix)
    c_lambdaRange = _c_arr(c_double, 2, lambdaRange)
    if roundNegativeCapacity:
        c_roundNegativeCapacity = c_int(1)
    else:
        c_roundNegativeCapacity = c_int(0)

    return {
        "numNodes": c_numNodes,
        "numArcs": c_numArcs,
        "source": c_source,
        "sink": c_sink,
        "arcMatrix": c_arcMatrix,
        "lambdaRange": c_lambdaRange,
        "roundNegativeCapacity": c_roundNegativeCapacity,
    }


def _create_c_output():
    c_numBreakpoints = c_int(0)
    c_cuts = POINTER(c_int)()
    c_breakpoints = POINTER(c_double)()
    c_stats = _c_arr(c_int, 5, (0,) * 5)
    c_times = _c_arr(c_double, 3, (0.0,) * 3)

    return {
        "numBreakpoints": c_numBreakpoints,
        "cuts": c_cuts,
        "breakpoints": c_breakpoints,
        "stats": c_stats,
        "times": c_times,
    }


def _solve(c_input, c_output):

    hpf_solve = libhpf.hpf_solve
    hpf_solve.argtypes = [
        c_int,
        c_int,
        c_int,
        c_int,
        POINTER(c_double),
        c_double * 2,
        c_int,
        POINTER(c_int),
        POINTER(POINTER(c_int)),
        POINTER(POINTER(c_double)),
        c_int * 5,
        c_double * 3,
    ]

    hpf_solve(
        c_input["numNodes"],
        c_input["numArcs"],
        c_input["source"],
        c_input["sink"],
        cast(byref(c_input["arcMatrix"]), POINTER(c_double)),
        c_input["lambdaRange"],
        c_input["roundNegativeCapacity"],
        byref(c_output["numBreakpoints"]),
        byref(c_output["cuts"]),
        byref(c_output["breakpoints"]),
        c_output["stats"],
        c_output["times"],
    )


def _cleanup(c_output):
    libhpf.libfree(c_output["breakpoints"])
    libhpf.libfree(c_output["cuts"])


def _check_multipliers_sink_adjacent_negative(G, sink, mult_cap):
    for u in G.predecessors(sink):
        if G[u][sink][mult_cap] > 0:
            raise ValueError(
                "Sink adjacent arcs should have non-positive multipliers. Arc (%s, %s = sink) has a multiplier of %f. Please reverse graph."
                % (u, sink, G[u][sink][mult_cap])
            )


def _check_multipliers_source_adjacent_positive(G, source, mult_cap):
    for v in G.successors(source):
        if G[source][v][mult_cap] < 0:
            raise ValueError(
                "Source adjacent arcs should have non-negative multipliers. Arc (%s = source, %s) has a multiplier of %f. Please reverse graph."
                % (source, v, G[source][v][mult_cap])
            )


def _read_output(c_output, nodeNames):
    numBreakpoints = c_output["numBreakpoints"].value
    breakpoints = [c_output["breakpoints"][i] for i in range(numBreakpoints)]

    cuts = {}
    for i, node in enumerate(nodeNames):
        cuts[node] = [
            c_output["cuts"][len(nodeNames) * j + i] for j in range(numBreakpoints)
        ]

    info = {
        "numArcScans": c_output["stats"][0],
        "numMergers": c_output["stats"][1],
        "numPushes": c_output["stats"][2],
        "numRelabels": c_output["stats"][3],
        "numGap": c_output["stats"][4],
        "readDataTime": c_output["times"][0],
        "intializationTime": c_output["times"][1],
        "solveTime": c_output["times"][2],
    }

    return breakpoints, cuts, info


def hpf(
    G,
    source,
    sink,
    const_cap,
    mult_cap=None,
    lambdaRange=None,
    roundNegativeCapacity=False,
):

    if mult_cap:
        parametric = True
        _check_multipliers_sink_adjacent_negative(G, sink, mult_cap)
        _check_multipliers_source_adjacent_positive(G, source, mult_cap)
    else:
        parametric = False
        lambdaRange = [0.0, 0.0]

    nodeNames, nodeDict, arcMatrix = _get_arcmatrix(
        G, const_cap, mult_cap, source, sink
    )

    c_input = _create_c_input(
        G, nodeDict, source, sink, arcMatrix, lambdaRange, roundNegativeCapacity
    )
    c_output = _create_c_output()

    _solve(c_input, c_output)

    breakpoints, cuts, info = _read_output(c_output, nodeNames)

    _cleanup(c_output)

    if not parametric:
        breakpoints = [None]

    return breakpoints, cuts, info
