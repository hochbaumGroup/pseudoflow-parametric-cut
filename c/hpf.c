/*************************************************************************
 * Hochbaum's Pseudo-flow (HPF) Algorithm for Parametric Minimimum Cut   *
 * ***********************************************************************
 * The HPF algorithm for finding Minimum-cut in a graph is described in: *
 * [1] D.S. Hochbaum, "The Pseudoflow algorithm: A new algorithm for the *
 * maximum flow problem", Operations Research, 58(4):992-1009,2008.      *
 *                                                                       *
 * The algorithm was found to be fast in theory (see the above paper)    *
 * and in practice (see:                                                 *
 * [2] D.S. Hochbaum and B. Chandran, "A Computational Study of the      *
 * Pseudoflow and Push-relabel Algorithms for the Maximum Flow Problem,  *
 * Operations Research, 57(2):358-376, 2009.                             *
 *                                                                       *
 * and                                                                   *
 *                                                                       *
 * [3] B. Fishbain, D.S. Hochbaum, S. Mueller, "Competitive Analysis of  *
 * Minimum-Cut Maximum Flow Algorithms in Vision Problems,               *
 * arXiv:1007.4531v2 [cs.CV]                                             *
 *                                                                       *
 * The algorithm solves a parametric s-t minimum cut problem. The		 *
 * algorithm finds all breakpoints for which the source set of the		 *
 * minimum cut changes as a function of lambda in the range				 *
 * [lower bound, upper bound] by recursively concluding that the interval  *
 * contains 0, 1, or more breakpoints. If the interval contains more than*
 * 1 breakpoint, then the interval is split into two interval, each of   *
 * which contains at least one breakpoint.								 *
 *                                                                       *
 * Parametric cut/flow problems allow for a linear function with input   *
 * lambda on source or sink adjacent arcs. Arcs that are adjacent to	 *
 * source should be non-decreasing in lambda and sink adjacent arcs		 *
 * should be non-increasing in lambda. The algorithm is able to deal with*
 * the reverse configuration (non-increasing on source adjacent arcs and *
 * non-decreasing on sink adjacent arcs) by flipping source and sink and *
 * reversing the direction of the arcs.									 *
 *                                                                       *
 * Usage:																 *
 * 1. Compile hpf.c with a C-compiler (e.g. gcc)						 *
 * 2. To execute within bash environment:								 *
 *	 <name compiled hpf executable> <path input file> <path output file> *
 *                                                                       *
 * INPUT FILE                                                            *
 * **********                                                            *
 * The input file is assumed to be in a modified DIMACS format:	         *
 * c <comment lines>													 *
 * p <# nodes> <# arcs> <lower bound> <upper bound> <round if negative>  *
 * n <source node> s													 *
 * n <sink node> t														 *
 * a <from-node> <to-node> <constant capacity> <lambda multiplier>		 *
 * where the following conditions are satisfied:						 *
 * - Nodes are labeled 0 .. <# nodes> - 1								 *
 * - <lambda multiplier> is non-negative if <from-node> == <source node> *
 *		and <to-node> != <sink-node>									 *
 * - <lambda multiplier> is non-positive if <from-node> != <source node> *
 *		and <to-node> == <sink-node>									 *
 * - <lambda multiplier> is zero if <from-node> != <source node>		 *
 *		and <to-node> != <sink-node>									 *
 * - <lambda multiplier> can take any value if							 *
 *		<from-node> != <source node> and <to-node> != <sink-node>		 *
 * - <round if negative> takes value 1 if the any negative capacity arc  *
 *		should be rounded to 0, and it takes value 0 otherwise			 *
 *                                                                       *
 * OUTPUT FILE                                                           *
 * ***********                                                           *
 * The solver will generate the following output file:					 *
 * t <time (in sec) read data> <time (in sec) initialize> <time			 *
		(in sec) solve>													 *
 * s <# arc scans> <# mergers> <# pushes> <# relabels > <# gap >		 *
 * p <number of lambda intervals = k>									 *
 * l <lambda upperbound interval 1> ... <lambda upperbound interval k>   *
 * n <node-id> <sourceset indicator intval 1 > .. <indicator intval k>   *
 *                                                                       *
 * Set-up                                                                *
 * ******                                                                *
 * Uncompress the MatlabHPF.zip file into the Matlab's working directory *
 * The zip file contains the following files:                            *
 * hpf.c - source code                                                   *
 * hpf.m - Matlab's help file                                            *
 * hpf.mexmaci - The compiled code for Mac OS 10.0.5 (Intel)/ Matlab     *
 *               7.6.0.324 (R2008a).                                     *
 * hpf.mexw32  - The compiled code for Windows 7 / Matlab 7.11.0.584     *
 *               (R2010b).                                               *
 * demo_general - Short Matlab code that generates small network and     *
 *                computes the minimum flow                              *
 * demo_vision - Short Matlab code that loads a Multiview reconstruction *
 *               vision problem (see: [3]) and computes its minimum cut. *
 * gargoyle-smal.mat - The vision problem.                               *
 *                                                                       *
 * When using this code, please cite:                                    *
 * References [1], [2] and [3] above and:                                *
 * Q. Spaen, B. Fishbain and D.S. Hochbaum, "Hochbaum's Pseudo-flow C	 *
 * Implementation", http://riot.ieor.berkeley.edu/riot/Applications/     *
 * Pseudoflow/maxflow.html                                               *
 *************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "../core/libhpf.h"

static void readData(char *filename, int* numNodes, int* numArcs, int *source, int *sink, double ** arcMatrixPointer, double lambdaRange[2], int * roundNegativeCapacity)
/*************************************************************************
readData
*************************************************************************/
{
	char buffer[32768];

	/* define parameters */
	int arcCount = 0;
	int numRemovedArcs = 0;
	int currentNode;
	int isSourceAssigned = 0;
	int isSinkAssigned = 0;
	char sourceSinkIndicator;
	int from;
	int to;
	double constantCapacity;
	double multiplierCapacity;
	double * arcMatrix;
	double * arcMatrixNew;

	int i;

	// open input file
	FILE* f = fopen(filename, "r");
	if (f == NULL)
	{
		printf("I/O error while opening input file %s", filename);
		exit(0);
	}

	/* Read lines of input file */
	while (1)
	{
		if (fgets(buffer, sizeof buffer, f) != NULL)
		{
			switch (*buffer)
			{
			case 'p': /* initialize problem */
				sscanf(buffer, "p %d %d %lf %lf %d\n", numNodes, numArcs, &lambdaRange[0], &lambdaRange[1], roundNegativeCapacity);

				if ((arcMatrix = (double *)malloc(*numArcs * 4 * sizeof(double))) == NULL)
				{
					printf("Could not allocate memory.\n");
					exit(0);
				}

				break;
			case 'n':
				sscanf(buffer, "n %i %c\n", &currentNode, &sourceSinkIndicator);
				if (sourceSinkIndicator == 's')
				{
					/* check if source is valid */
					if (currentNode >= *numNodes || currentNode < 0)
					{
						printf("Nodes are labeled from 0 to <number of nodes>  - 1\n");
						exit(0);
					}
					/* check if source is assigned */
					if (isSourceAssigned)
					{
						printf("Source is already defined\n");
						exit(0);
					}
					else
					{
						*source = currentNode;
						isSourceAssigned = 1;
					}
				}
				else if (sourceSinkIndicator == 't')
				{
					/* check if sink is valid */
					if (currentNode >= *numNodes || currentNode < 0)
					{
						printf("Nodes are labeled from 0 to <number of nodes>  - 1\n");
						exit(0);
					}
					/* check if sink is assigned */
					if (isSinkAssigned)
					{
						printf("Sink is already defined\n");
						exit(0);
					}
					else
					{
						*sink = currentNode;
						isSinkAssigned = 1;
					}
				}
				else
				{
					printf("Node type: %c is unknown\n", sourceSinkIndicator);
					exit(0);
				}

				break;
			case 'a':
				sscanf(buffer, "a %d %d %lf %lf\n", &from, &to, &constantCapacity, &multiplierCapacity);

				if (isSinkAssigned == 0 || isSourceAssigned == 0)
				{
					printf("Source and sink need to be defined before arcs are defined.\n");
					exit(0);
				}

				/* assign arc */
				if (from < 0 || to < 0 || from >= *numNodes || to >= *numNodes )
				{
					printf("Nodes are labeled from 0 to <number of nodes>  - 1\n");
					exit(0);
				}
				else if (from == to)
				{
					printf("Node %u has a self loop which is not allowed\n", from);
					exit(0);
				}
				else if (multiplierCapacity > 0 && from != *source)
				{
					printf("Only source adjacent arcs can have a strictly positive capacity multiplier\n");
					exit(0);
				}
				else if (multiplierCapacity < 0 && to != *sink)
				{
					printf("Only sink adjacent arcs can have a strictly negative capacity multiplier\n");
					exit(0);
				}
				else if (arcCount >= *numArcs)
				{
					printf("Incorrect number of arcs specified\n");
					exit(0);
				}
				else if (to==*source || from==*sink)
				{
				  numRemovedArcs++;
				  continue;
				}

				arcMatrix[ arcCount * 4 + 0] = (double) from;
				arcMatrix[ arcCount * 4 + 1] = (double) to;
				arcMatrix[ arcCount * 4 + 2] = constantCapacity;
				arcMatrix[ arcCount * 4 + 3] = multiplierCapacity;

				++arcCount;
			}
		}
		else if (feof(f))
		{
			break;
		}
		else
		{
			printf("I/O error while reading %s\n", filename);
			exit(0);
		}
	}

	// close file
	fclose(f);

	if (numRemovedArcs>0)
	{
	  *numArcs-= numRemovedArcs;
	  if ((arcMatrixNew = (double *)malloc(*numArcs * sizeof(double))) == NULL)
	    {
	      printf("Could not allocate memory.\n");
	      exit(0);
	    }
	  for(i = 0;i < *numArcs * 4;++i)
		{
	   	arcMatrixNew[i] = arcMatrix[i];
	  }
	  free(arcMatrix);
	  arcMatrix = arcMatrixNew;
	}
	/* check if correct number of arcs has been specified */
	if (arcCount != *numArcs)
	{
		printf("Incorrect number of arcs specified\n");
		exit(0);
	}
	else if (isSourceAssigned == 0)
	{
		printf("Source is not assigned\n");
		exit(0);
	}
	else if (isSinkAssigned == 0)
	{
		printf("Sink is not assigned\n");
		exit(0);
	}
	else if (*source == *sink)
	{
		printf("The source node and sink node need to be distinct\n");
		exit(0);
	}

	*arcMatrixPointer = arcMatrix;
}

static void writeOutput (char *filename, int numBreakpoints, int numNodes, double* breakpoints, int* cuts, int* stats, double* times)
{
/*************************************************************************
writeOutput
*************************************************************************/
	int i, j;

	// open outputFile
	FILE* f = fopen(filename,"w");
	if (f == NULL)
	{
		printf("I/O error while opening output file %s", filename);
		exit(0);
	}

	/* print times */
	fprintf(f, "t %.3lf %.3lf %.3lf\n", times[0], times[1], times[2]);
	/* print stats */
	fprintf(f, "s %d %d %d %d %d\n", stats[0], stats[1], stats[2], stats[3], stats[4]);

	fprintf(f, "p %d\n", numBreakpoints);

	/* print lambda values */
	fprintf(f, "l ");
	for (i = 0; i < numBreakpoints; i++)
	{
		fprintf(f, "%.12lf", breakpoints[i]);
		if (i < numBreakpoints - 1)
		{
			fprintf(f, " ");
		}
		else
		{
			fprintf(f, "\n");
		}
	}

	/* print values nodes*/
	for (i = 0; i < numNodes; i++)
	{
		fprintf(f, "n %d ",i);
		for (j = 0; j < numBreakpoints; j++)
		{
			fprintf(f, "%d", cuts[i + j * numNodes]);
			if (j < numBreakpoints - 1)
			{
				fprintf(f, " ");
			}
			else
			{
				fprintf(f, "\n");
			}
		}
	}

	// close output file
	fclose(f);
}


int main(int argc, char **argv)
/*************************************************************************
main - Main function
*************************************************************************/
{
	// check number of input arguments
	if (argc != 3)
	{
		printf("Incorrect number of input arguments. Call hpf.exe inputFile outputFile\n");
		exit(0);
	}

	// prepare input solver
	int numNodes;
	int numArcs;
	int source;
	int sink;
	double* arcMatrix;
	double lambdaRange[2];
	int roundNegativeCapacity;

	readData(argv[1], &numNodes, &numArcs, &source, &sink, &arcMatrix, lambdaRange, &roundNegativeCapacity);

	printf("NumNodes: %d\n", numNodes);
	printf("NumArcs: %d\n", numArcs);
	printf("Lambda Range: [%lf, %lf]\n", lambdaRange[0], lambdaRange[1]);
	printf("Round if negative: %d\n", roundNegativeCapacity);
	printf("Arc matrix:\n");
	for (int i = 0; i < numArcs; ++i)
	{
		printf("Row %d: [%.2lf, %.2lf, %.2lf, %.2lf]\n", i, arcMatrix[i * 4 + 0 ], arcMatrix[i * 4 + 1 ], arcMatrix[i * 4 + 2 ], arcMatrix[i * 4 + 3 ]);
	}


	// prepare output solver
	int numBreakpoints;
	int *cuts;
	double *breakpoints;
	int stats[5];
	double times[3];

	hpf_solve(numNodes, numArcs, source, sink, arcMatrix, lambdaRange, roundNegativeCapacity, &numBreakpoints, &cuts, &breakpoints, stats, times );

	printf("Stats: [%d, %d, %d, %d, %d]\n", stats[0],stats[1],stats[2],stats[3],stats[4]);
	printf("times: [%lf, %lf, %lf]\n", times[0],times[1],times[2]);
	printf("Num breakpoints: %d\n", numBreakpoints);
	printf("breakpoints:\n");
	for (int i = 0; i < numBreakpoints; ++i)
	{
		printf("Breakpoint: %lf\n", breakpoints[i]);
		for (int j = 0; j < numNodes; j++)
		{
			printf("Cut indicator: %d\n", cuts[i * (int) numNodes+ j ]);
		}
	}

	writeOutput(argv[2], numBreakpoints, numNodes, breakpoints, cuts, stats, times);

	return 1;
}
