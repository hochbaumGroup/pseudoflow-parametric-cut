import unittest
from pseudoflow.python import hpf
from ctypes import c_int
from networkx import DiGraph


def get_vals_array(array, size):
    return [array[i] for i in range(size)]


class TestPythonHpf(unittest.TestCase):

    def test_c_arr(self):
        arr = hpf._c_arr(c_int, 10, range(10))
        TenIntArr = c_int * 10
        self.assertEqual(get_vals_array(arr, 10),
                         get_vals_array(TenIntArr(0, 1, 2, 3, 4, 5,
                                                  6, 7, 8, 9),
                                        10)
                         )

    def test_get_arc_matrix(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])

        G[0][1]["const"] = 1
        G[1][2]["const"] = 9.0

        G[0][1]["mult"] = 5
        G[1][2]["mult"] = -3

        linArcMat = [0., 1., 1., 5., 1., 2., 9., -3.]
        nodeNames, nodeDict, arcMatOut = hpf._get_arcmatrix(G,
                                                            const_cap="const",
                                                            mult_cap="mult")
        self.assertEqual(nodeNames, [0, 1, 2])
        self.assertEqual(nodeDict, {0: 0, 1: 1, 2: 2})
        self.assertEqual(list(arcMatOut), linArcMat)

    def test_create_c_input(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])

        G[0][1]["const"] = 1
        G[1][2]["const"] = 9.0

        G[0][1]["mult"] = 5
        G[1][2]["mult"] = -3

        source = 0
        sink = 2
        roundNegativeCapacity = True
        lambdaRange = [0., 2.]

        linArcMat = [0., 1., 1., 5., 1., 2., 9., -3.]
        nodeDict = {0: 0, 1: 1, 2: 2}

        c_input = hpf._create_c_input(G, nodeDict, source, sink,
                                      linArcMat, lambdaRange,
                                      roundNegativeCapacity)

        self.assertEqual(c_input['numNodes'].value, 3)
        self.assertEqual(c_input['numArcs'].value, 2)
        self.assertEqual(c_input['source'].value, 0)
        self.assertEqual(c_input['sink'].value, 2)
        self.assertEqual(len(c_input['arcMatrix']), 2 * 4)
        self.assertEqual(get_vals_array(c_input['arcMatrix'], 8),
                         linArcMat)
        self.assertEqual(get_vals_array(c_input['lambdaRange'], 2),
                         lambdaRange)
        self.assertEqual(c_input['roundNegativeCapacity'].value, 1)

    def test_create_c_output(self):
        c_output = hpf._create_c_output()
        self.assertEqual(c_output['numBreakpoints'].value, 0)
        self.assertEqual(len(c_output['stats']), 5)
        self.assertEqual(get_vals_array(c_output['stats'], 5),
                         [0., 0., 0., 0., 0.])
        self.assertEqual(len(c_output['times']), 3)
        self.assertEqual(get_vals_array(c_output['times'], 3),
                         [0., 0., 0.])

    def test_hpf_solve_nonparametric(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])

        G[0][1]["const"] = 1
        G[1][2]["const"] = 9.0

        G[0][1]["mult"] = 5
        G[1][2]["mult"] = -3

        source = 0
        sink = 2
        roundNegativeCapacity = True
        lambdaRange = [0., 0.]

        nodeNames, nodeDict, linArcMat = hpf._get_arcmatrix(G, "const", "mult")

        c_input = hpf._create_c_input(G, nodeDict, source, sink,
                                      linArcMat, lambdaRange,
                                      roundNegativeCapacity)
        c_output = hpf._create_c_output()

        hpf._solve(c_input, c_output)

        self.assertEqual(c_output['numBreakpoints'].value, 1)
        self.assertEqual(get_vals_array(c_output['breakpoints'], 1),
                         [0.])
        self.assertEqual(get_vals_array(c_output['cuts'], 3),
                         [1, 0, 0])

    def test_hpf_solve_parametric(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])

        G[0][1]["const"] = 1
        G[1][2]["const"] = 9.0

        G[0][1]["mult"] = 5
        G[1][2]["mult"] = -3

        source = 0
        sink = 2
        roundNegativeCapacity = True
        lambdaRange = [0., 2.]

        nodeNames, nodeDict, linArcMat = hpf._get_arcmatrix(G, "const", "mult")

        c_input = hpf._create_c_input(G, nodeDict, source, sink,
                                      linArcMat, lambdaRange,
                                      roundNegativeCapacity)
        c_output = hpf._create_c_output()

        hpf._solve(c_input, c_output)

        self.assertEqual(c_output['numBreakpoints'].value, 2)
        self.assertEqual(get_vals_array(c_output['breakpoints'], 2),
                         [1., 2.])
        self.assertEqual(get_vals_array(c_output['cuts'], 6),
                         [1, 0, 0, 1, 1, 0])

    def test_read_output(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])

        G[0][1]["const"] = 1
        G[1][2]["const"] = 9.0

        G[0][1]["mult"] = 5
        G[1][2]["mult"] = -3

        source = 0
        sink = 2
        roundNegativeCapacity = True
        lambdaRange = [0., 2.]

        nodeNames, nodeDict, linArcMat = hpf._get_arcmatrix(G, "const", "mult")

        c_input = hpf._create_c_input(G, nodeDict, source, sink,
                                      linArcMat, lambdaRange,
                                      roundNegativeCapacity)
        c_output = hpf._create_c_output()

        hpf._solve(c_input, c_output)

        breakpoints, cuts, info = hpf._read_output(c_output, nodeNames)

        self.assertEqual(breakpoints, [1., 2.])
        self.assertEqual(cuts, {0: [1, 1], 1: [0, 1], 2: [0, 0]})

        # TODO: Implement test for info! Results are currently incorrect.

    def test_hpf(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])

        G[0][1]["const"] = 1
        G[1][2]["const"] = 9.0

        G[0][1]["mult"] = 5
        G[1][2]["mult"] = -3

        source = 0
        sink = 2
        roundNegativeCapacity = True
        lambdaRange = [0., 2.]

        breakpoints, cuts, info = hpf.hpf(
            G, source, sink, const_cap="const", mult_cap="mult",
            lambdaRange=lambdaRange, roundNegativeCapacity=False)

        self.assertEqual(breakpoints, [1., 2.])
        self.assertEqual(cuts, {0: [1, 1], 1: [0, 1], 2: [0, 0]})

        # TODO: Implement test for info! Results are currently incorrect.
