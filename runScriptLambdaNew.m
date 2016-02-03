% compile mex file for running c Library in Matlab
mex -largeArrayDims hpf.c

lowProblem = CutProblem( 8, 1, 8 , capLabels, cap, sourceWeights, sinkWeights, lambdaMultSource, lambdaMultSink, 0, 0, lambdaLow );
highProblem = CutProblem( 8, 1, 8 , capLabels, cap, sourceWeights, sinkWeights, lambdaMultSource, lambdaMultSink, 0, 0, lambdaHigh );

[ lambdas, cuts ] = hpfCompleteParametric( lowProblem, highProblem );