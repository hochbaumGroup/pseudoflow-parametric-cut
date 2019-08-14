# Hochbaum's Pseudoflow (HPF) Algorithm for (Linear) Fully Parametric Minimum Cut
This package provides a parametric implementation of pseudoflow for minimum cut on directed graphs. In the parametric minimum cut problem, the capacity of source-adjacent arcs is monotone non-decreasing in the parameter `lambda` whereas the capacity of sink-adjacent arcs is monotone non-increasing in `lambda`. This solver requires that the capacities of source and sink adjacent arcs are linear in `lambda`: `capacity = constant + multiplier * lambda`.

This fully parametric solver finds the optimal minimum cut for all `lambda` values in a given range. The solution for all lambda values is represented with `O(n)` intervals for the parameter lambda. In each interval, the optimal minimum cut remains the same.

A simple parametric minimum cut solver that provides the optimal minimum cut for a given list of arc capacities is available [here](https://riot.ieor.berkeley.edu/Applications/Pseudoflow/parametric.html), and a non-parametric maximum flow version of pseudoflow is available [here](https://riot.ieor.berkeley.edu/Applications/Pseudoflow/maxflow.html).

The package provides interfaces for Python, C, and Matlab.

This implementation uses a variant of the fully parametric HPF algorithm as described in:
>    DS Hochbaum (2008), The Pseudoflow algorithm: A new algorithm for the maximum flow problem. Operations Research, 58(4):992-1009.

This implementation does not use *free runs* nor does it use warm starts with informatiom from previous runs (see pg.15). This implementation should therefore **not be used** for comparison with the fully parametric HPF algorithm.

The package provides an option to round capacities that are negative for certain lambda values to zero. This option should **only** be used when each node has a source adjacent arc with capacity `max(0, a * lambda + b)` and a corresponding sink adjacent arc with capacity `max(0, -a * lambda - b)`. Otherwise, the intersection of the cut capacities is wrongly identified.


## Instructions for Python

Install the package with `pip`:

```bash
    pip install pseudoflow
```

#### Example
```python
import networkx as nx
import pseudoflow

G = nx.DiGraph()
G.add_edge(0, 1, const=1, mult=5)
G.add_edge(1, 2, const=9, mult=-3)


source = 0
sink = 2
lambda_range = [0., 2.]

breakpoints, cuts, info = pseudoflow.hpf(
    G,  # Networkx directed graph.
    source,  # Node id of the source node.
    sink,  # Node id of the sink node.
    const_cap="const",  # Edge attribute with the constant capacity.
    mult_cap="mult",  # Edge attribute with the lambda multiplier.
    lambdaRange=lambda_range,  # (lower, upper) bounds for the lambda parameter.
    roundNegativeCapacity=False  # True if negative arc capacities should be rounded to zero.
)

# breakpoints: list of upper bounds for the lambda intervals.
# cuts: A dictionary with for each node a list indicating whether
#       the node is in the source set of the minimum cut.
print(breakpoints)  # Output: [1., 2.]
print(cuts)  # Output: {0: [1, 1], 1: [0, 1], 2: [0, 0]}
```

## Instructions for C
Navigate to directory `src/pseudoflow/c`, and compile the `hpf` executable with `make`.

To execute the solver, use:
```bash
hpf input-file.txt output-file.txt
```

The input file should contain the graph structure and is assumed to have the following format:
```
    c <comment lines>
    p <# nodes> <# arcs> <lower bound> <upper bound> <round if negative>
    n <source node> s
    n <sink node> t
    a <from-node> <to-node> <constant capacity> <lambda multiplier>
```
where the `a` line is repeated for each arc. The file should satisfy the following conditions:
- Nodes are labeled `0 .. <# nodes> - 1`.
- `<lambda multiplier>` is non-negative if `<from-node> == <source node>` and `<to-node> != <sink-node>`.
- `<lambda multiplier>` is non-positive if `<from-node> != <source node>` and `<to-node> == <sink-node>`.
- `<lambda multiplier>` is zero if `<from-node> != <source node>` and `<to-node> != <sink-node>`.
- `<round if negative>` takes value 1 if the any negative capacity arc should be rounded to 0, and value 0 otherwise.

The solver will generate the following output file:
```
t <time (in sec) read data> <time (in sec) initialize> <time (in sec) solve>
s <# arc scans> <# mergers> <# pushes> <# relabels > <# gap >
p <number of lambda intervals = k>
l <lambda upperbound interval 1> ... <lambda upperbound interval k>
n <node-id> <sourceset indicator interval 1 > .. <indicator interval k>
```
The `n` line appears for each node. `<sourceset indicator interval 1 >` indicates whether the node is in the source set of the minimum cut for the first lambda interval.

See `src/pseudoflow/c/example` for an example.

## Instructions for Matlab

Copy the content of `src/pseudoflow/matlab` to your current directory.

From within Matlab, compile the mex extension with:
```matlab
    mex hpfMatlab.c
```

The solver is accessible via the `hpf` function with the following signature:
```matlab
    [cuts, lambdas, stats, times]  = hpf(arcmatrix, num_nodes, source, sink lambda_range, rounding);
```

#### Inputs:
* **arcmatrix**: Each row of the matrix has the following structure: `[from_node, to_node, constant capacity, lambda multiplier]`
* **num_nodes**: Number of nodes in the graph
* **source_node**: The numeric label of the source node
* **sink_node**: The numeric label of the sink node
* **lambda_range**: [lower bound, upper bound] for the lambda parameter.
* **rounding**: Set to 1 if negative arc capacities should be rounded to zero, and 0 otherwise.

#### Outputs:
* **cuts**: n x k matrix where `A(i,j)` is 1 if node `i` is in the source set for lambda interval `j`, and 0 otherwise.
* **lambdas**: 1 x k matrix where `L(j)` is the upper bound of the lambda interval `j`.
