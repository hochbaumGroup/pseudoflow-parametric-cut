void hpf_solve(int numNodes, int numArcs, int source, int sink, double * arcMatrix, double lambdaRange[2], int roundNegativeCapacityIn, int * numBreakpoints, int ** cuts, double ** breakpoints, int stats[5], double times[3] );

void libfree(void * p);
