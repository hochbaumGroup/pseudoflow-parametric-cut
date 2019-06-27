function [cuts, lambdas] = hpf( arcmatrix, num_nodes, source_node, sink_node, lambda_range, rounding );
%   Hochbaum's Pseudo-flow (HPF) Algorithm Matlab implementation
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The HPF algorithm for finding Minimum-cut in a graph is described in:
%   [1] D.S. Hochbaum, "The Pseudoflow algorithm: A new algorithm for the
%   maximum flow problem", Operations Research, 58(4):992-1009,2008.
%
%   The algorithm was found to be fast in theory (see the above paper)
%   and in practice (see:
%   [2] D.S. Hochbaum and B. Chandran, "A Computational Study of the
%   Pseudoflow and Push-relabel Algorithms for the Maximum Flow Problem,
%   Operations Research, 57(2):358-376, 2009.
%
%   and
%
%   [3] B. Fishbain, D.S. Hochbaum, S. Mueller, "Competitive Analysis of
%   Minimum-Cut Maximum Flow Algorithms in Vision Problems,
%   arXiv:1007.4531v2 [cs.CV]
%
%   Usage: Within Matlab environment:
%   [cuts, lambdas, stats, times]  = hpf(arcmatrix, num_nodes, source, sink lambda_range, rounding);
%
%   INPUTS
%
%   arcmatrix - Each row of the matrix is as follows:
%               [from_node, to_node, constant capacity, lambda multiplier]
%   num_nodes - Number of nodes in the graph
%   source_node - The numeric label of the source node
%   sink_node   - The numeric label of the sink node
%   lambda_range - [lower bound, upper bound] for lambda values.
%   rounding - 1 if negative arc capacities should be rounded to zero, and 0 otherwise.
%
%
%   OUTPUTS
%
%   cuts - n x k matrix where A(i,j) is 1 if node i is in the source cet for lambda interval j.
%   lambdas - 1 x k matrix where L(j) is the upper bound of the j-th lambda interval.

[cuts, lambdas, stats, time] = hpfMatlab( arcmatrix, num_nodes, source_node, sink_node, lambda_range, rounding );
