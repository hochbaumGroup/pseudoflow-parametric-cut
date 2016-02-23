% -------------------------------------------------------------------------
% Example Parametric cut
% Author: Quico Spaen
% -------------------------------------------------------------------------
% Solve a parametric minimum cut problem from 1 (source) to 8 (sink) for 
% lambda in range [0,2]. It correctly identifies breakpoints: [1, 4/3] 
% (2 is not a break point) and displays the cut for ranges in [0,1], 
% (1,4/3], (4/3, 2]. Note that the lambdas vector gives the upper bound 
% for each lambda interval.
%
% The arcs in the graph are as follows:
% (1,2) = 7
% (1,3) = 7
% (2,4) = 5
% (2,5) = 2
% (3,2) = 3
% (3,4) = 4
% (3,7) = 5
% (4,5) = 2
% (4,6) = 3
% (4,7) = 4
% (5,8) = 10 - 1.5 lambda
% (6,5) = 4
% (6,8) = 2 - lambda
% (7,6) = 2
% (7,8) = 3

% compile mex file for running c Library in Matlab
mex -largeArrayDims hpf.c

% load data
load('example.mat')

% instantiated cut instances. One with lambda = 0 (lower bound) and one
% with lambda = 2 (upper bound)
lowProblem = CutProblem( 8, 1, 8 , capLabels, cap, sourceWeights, sinkWeights, lambdaMultSource, lambdaMultSink, 0, 0, lambdaLow );
highProblem = CutProblem( 8, 1, 8 , capLabels, cap, sourceWeights, sinkWeights, lambdaMultSource, lambdaMultSink, 0, 0, lambdaHigh );

[ lambdas, cuts ] = hpfCompleteParametric( lowProblem, highProblem );