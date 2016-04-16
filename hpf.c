/*************************************************************************
 * Hochbaum's Pseudo-flow (HPF) Algorithm Matlab implementation          *
 * ************************************************************          *
 * The HPF algorithm for finding Minimum-cut in a graph is described in: *                                        *
 * [1] D.S. Hochbaum, "The Pseudoflow algorithm: A new algorithm for the *
 * maximum flow problem", Operations Research, 58(4):992-1009,2008.      *
 *                                                                       *
 * The algorithm was found to be fast in theory (see the above paper)    *
 * and in practice (see:                                                 *
 * [2] D.S. Hochbaum and B. Chandran, "A Computational Study of the      *
 * Pseudoflow and Push-relabel Algorithms for the Maximum Flow Problem,  *
 * Operations Research, 57(2):358-376, 2009.                             *
 *
 * and                                                                   *
 *                                                                       *
 * [3] B. Fishbain, D.S. Hochbaum, S. Mueller, "Competitive Analysis of  *
 * Minimum-Cut Maximum Flow Algorithms in Vision Problems,               *
 * arXiv:1007.4531v2 [cs.CV]                                             *
 *                                                                       *
 * Usage: Within Matlab environment:                                     *
 * [value,cut] = hpf(sim_mat,source,sink);                               *
 *                                                                       *
 * INPUTS                                                                *
 * ******                                                                *
 * sim_mat - similarity matrix - a_{i,j} is the capacity of the arc (i,j)*
 *           a_{i,j} are non-negatives; the self-similarities (diagonal  *
 *           values) are zero.                                           *
 *           the sim_mat should be sparse (see Matlab's help)            *
 * source - The numeric label of the source node                         *
 * sink   - The numeric label of the sink node                           *
 *                                                                       *
 * OUTPUTS                                                               *
 * *******                                                               *
 * value - the capacity of the cut                                       *
 * cut   - the source set (see [1]), where x_i = 1, if i \in S ; 0 o/w   *
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
 * B. Fishbain and D.S. Hochbaum, "Hochbaum's Pseudo-flow Matlab         *
 * implementation", http://riot.ieor.berkeley.edu/riot/Applications/     *
 * Pseudoflow/maxflow.html                                               *
 *************************************************************************/

//#include "matrix.h"
#include "stdio.h"
#include "assert.h"
//#include <sys/time.h>
//#include <sys/resource.h>
#include "stdlib.h"
#include "time.h"
//#include <unistd.h>

/*************************************************************************
Definitions
*************************************************************************/
#define  MAX_LEVELS  300
#define VERSION 3.3

typedef unsigned int uint;
typedef long int lint;
typedef long long int llint;
typedef unsigned long long int ullint;

typedef struct cutProblem
{
	double lambdaValue;
	/* Arc list */
	/* Nod list */
	double cutValue;
	double cutMultiplier;
	double cutConstant;
	/* source set */
	/* sink set */
	double sourceSinkArcCapacity;
	double sourceSinkArcConstant;
	double sourceSinkArcMultiplier;
} CutProblem;

typedef struct
{
	uint label;
	uint numAdjacent;
} NodeSuper;

typedef struct
{
	NodeSuper *from;
	NodeSuper *to;
	double capacity;
	double constant;
	double multiplier;
} ArcSuper;

typedef struct arc
	{
		struct node *from;
		struct node *to;
		double flow;
		double capacity;
		uint direction;
	} Arc;

typedef struct node
	{
		uint visited;
		uint numAdjacent;
		uint number;
		uint label;
		double excess;
		struct node *parent;
		struct node *childList;
		struct node *nextScan;
		uint numOutOfTree;
		Arc **outOfTree;
		uint nextArc;
		Arc *arcToParent;
		struct node *next;
	} Node;

typedef struct root
	{
		Node *start;
		Node *end;
	} Root;

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

/*************************************************************************
Global variables
*************************************************************************/
static uint numNodes = NULL;
static uint numArcs = NULL;
static uint source;
static uint sink;
static uint highestStrongLabel = 1;

static double numArcScans = 0;
static double numPushes = 0;
static double numMergers = 0;
static double numRelabels = 0;
static double numGaps = 0;

static Node *nodesList = NULL;
static Root *strongRoots = NULL;
static uint *labelCount = NULL;
static Arc *arcList = NULL;
static NodeSuper *nodeListSuper = NULL;
static ArcSuper *arcListSuper = NULL;
static uint lowestPositiveExcessNode = 0;

static uint useParametricCut = 0;
static uint roundNegativeWeight = NULL;

double dabs(double value)
{
	if (value >= 0)
		return value;
	else return -value;
}

int isFlow(double flow)
/*************************************************************************
isFlow: We set a threshhold. If the flow value is below the threshhold, we
take it as no flow. Otherwise we take it as a flow
*************************************************************************/
{
	if (flow > 0)
		return 1;
	else return 0;
}

int isExcess(double excess)
/*************************************************************************
isExcess: We set a threshhold. If the absolute value of the excess is within
the threshold, then we take it as nothing. Otherwise we will return the sign
of the excess (deficit if negative).
*************************************************************************/
{
	if (excess < 0)
		return -1;
	else if (excess > 0)
		return 1;
	else return 0;
}

static void createOutOfTree (Node *nd)
{
/*************************************************************************
createOutOfTree
*************************************************************************/
	if (nd->numAdjacent)
	{
		if ((nd->outOfTree = (Arc **) malloc (nd->numAdjacent * sizeof (Arc *))) == NULL)
		{
			printf("Out of memory\n");
			exit(0);
		}
	}
}

static void initializeArc (Arc *ac)
{
/*************************************************************************
initializeArc
*************************************************************************/
	ac->from = NULL;
	ac->to = NULL;
	ac->capacity = 0;
	ac->flow = 0;
	ac->direction = 1;
}

static void initializeArcSuper(ArcSuper *ac)
{
	/*************************************************************************
	initializeArcSuper
	*************************************************************************/
	ac->from = NULL;
	ac->to = NULL;
	ac->capacity = 0;
	ac->constant = 0;
	ac->multiplier = 0;
}

static void liftAll (Node *rootNode)
{
/*************************************************************************
liftAll
*************************************************************************/
	Node *temp, *current=rootNode;

	current->nextScan = current->childList;

	-- labelCount[current->label];
	current->label = numNodes;

	for ( ; (current); current = current->parent)
	{
		while (current->nextScan)
		{
			temp = current->nextScan;
			current->nextScan = current->nextScan->next;
			current = temp;
			current->nextScan = current->childList;

			-- labelCount[current->label];
			current->label = numNodes;
		}
	}
}
static void addOutOfTreeNode (Node *n, Arc *out)
{
/*************************************************************************
addOutOfTreeNode
*************************************************************************/
	n->outOfTree[n->numOutOfTree] = out;
	++ n->numOutOfTree;
}

static void addToStrongBucket (Node *newRoot, Root *rootBucket)
{
/*************************************************************************
addToStrongBucket
*************************************************************************/
	if (rootBucket->start)
	{
		rootBucket->end->next = newRoot;
		rootBucket->end = newRoot;
		newRoot->next = NULL;
	}
	else
	{
		rootBucket->start = newRoot;
		rootBucket->end = newRoot;
		newRoot->next = NULL;
	}
}

static __inline int addRelationship (Node *newParent, Node *child)
{
/*************************************************************************
addRelationship
*************************************************************************/
	child->parent = newParent;
	child->next = newParent->childList;
	newParent->childList = child;

	return 0;
}

static __inline void breakRelationship (Node *oldParent, Node *child)
{
/*************************************************************************
breakRelationship
*************************************************************************/
	Node *current;

	child->parent = NULL;

	if (oldParent->childList == child)
	{
		oldParent->childList = child->next;
		child->next = NULL;
		return;
	}

	for (current = oldParent->childList; (current->next != child); current = current->next);

	current->next = child->next;
	child->next = NULL;
}

static void merge (Node *parent, Node *child, Arc *newArc)
{
/*************************************************************************
merge
*************************************************************************/
	Arc *oldArc;
	Node *current = child, *oldParent, *newParent = parent;

	++ numMergers;

	while (current->parent)
	{
		oldArc = current->arcToParent;
		current->arcToParent = newArc;
		oldParent = current->parent;
		breakRelationship (oldParent, current);
		addRelationship (newParent, current);
		newParent = current;
		current = oldParent;
		newArc = oldArc;
		newArc->direction = 1 - newArc->direction;
	}

	current->arcToParent = newArc;
	addRelationship (newParent, current);
}


static __inline void pushUpward (Arc *currentArc, Node *child, Node *parent, const double resCap)
{
/*************************************************************************
pushUpward
*************************************************************************/
	++ numPushes;

	if (isExcess(resCap-child->excess) >= 0)//(/*(int)*/resCap >= child->excess)
	{
		parent->excess += child->excess;
		currentArc->flow += child->excess;
		child->excess = 0;
		return;
	}

	currentArc->direction = 0;
	parent->excess += resCap;
	child->excess -= resCap;
	currentArc->flow = currentArc->capacity;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, &strongRoots[child->label]);
}


static __inline void pushDownward (Arc *currentArc, Node *child, Node *parent, double flow)
{
/*************************************************************************
pushDownward
*************************************************************************/
	++ numPushes;

	if (isExcess(flow - child->excess) >= 0)//(/*(int)*/flow >= child->excess)
	{
		parent->excess += child->excess;
		currentArc->flow -= child->excess;
		child->excess = 0;
		return;
	}

	currentArc->direction = 1;
	child->excess -= flow;
	parent->excess += flow;
	currentArc->flow = 0;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, &strongRoots[child->label]);
}

static void pushExcess (Node *strongRoot)
{
/*************************************************************************
pushExcess
*************************************************************************/
	Node *current, *parent;
	Arc *arcToParent;
	/*int*/double prevEx=1;

	for (current = strongRoot; (isExcess(current->excess) && current->parent); current = parent)
	{
		parent = current->parent;
		prevEx = parent->excess;

		arcToParent = current->arcToParent;

		if (arcToParent->direction)
		{
			pushUpward (arcToParent, current, parent, (arcToParent->capacity - arcToParent->flow));
		}
		else
		{
			pushDownward (arcToParent, current, parent, arcToParent->flow);
		}
	}

	if ((isExcess(current->excess) > 0) && (isExcess(prevEx) <= 0))
	{
		addToStrongBucket (current, &strongRoots[current->label]);
	}
}


static Arc * findWeakNode (Node *strongNode, Node **weakNode)
{
/*************************************************************************
findWeakNode
*************************************************************************/
	uint i, size;
	Arc *out;

	size = strongNode->numOutOfTree;

	for (i=strongNode->nextArc; i<size; ++i)
	{
		++ numArcScans;
		if (strongNode->outOfTree[i]->to->label == (highestStrongLabel-1))
		{
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->to;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		} else if (strongNode->outOfTree[i]->from->label == (highestStrongLabel-1)) {
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->from;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		}
	}

	strongNode->nextArc = strongNode->numOutOfTree;

	return NULL;
}


static void checkChildren (Node *curNode)
{
/*************************************************************************
checkChildren
*************************************************************************/
	for ( ; (curNode->nextScan); curNode->nextScan = curNode->nextScan->next)
	{
		if (curNode->nextScan->label == curNode->label)
		{
			return;
		}

	}

	-- labelCount[curNode->label];
	++	curNode->label;
	++ labelCount[curNode->label];

	++numRelabels;

	curNode->nextArc = 0;
}


static void simpleInitialization (void)
{
/*************************************************************************
simpleInitialization
*************************************************************************/
	uint i, size;
	Arc *tempArc;

	size = nodesList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) // Saturating source adjacent nodes
	{
		tempArc = nodesList[source-1].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->to->excess += tempArc->capacity;
	}

	size = nodesList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i) // Pushing maximum flow on sink adjacent nodes
	{
		tempArc = nodesList[sink-1].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->from->excess -= tempArc->capacity;
	}

	nodesList[source-1].excess = 0; // zeroing source excess
	nodesList[sink-1].excess = 0;	// zeroing sink excess

	for (i=0; i<numNodes; ++i)
	{
		if (isExcess(nodesList[i].excess) > 0)
		{
		    nodesList[i].label = 1;
			++ labelCount[1];

			addToStrongBucket (&nodesList[i], &strongRoots[1]);
		}
	}

	nodesList[source-1].label = numNodes;	// Set the source label to n
	nodesList[sink-1].label = 0;			// set the sink label to 0
	labelCount[0] = (numNodes - 2) - labelCount[1];
}


static Node* getHighestStrongRoot (void)
{
/*************************************************************************
getHighestStrongRoot
*************************************************************************/
	uint i;
	Node *strongRoot;

	for (i=highestStrongLabel; i>0; --i)
	{
		if (strongRoots[i].start)
		{
			highestStrongLabel = i;
			if (labelCount[i-1])
			{
				strongRoot = strongRoots[i].start;
				strongRoots[i].start = strongRoot->next;
				strongRoot->next = NULL;
				return strongRoot;
			}

			while (strongRoots[i].start)
			{
				++ numGaps;

				strongRoot = strongRoots[i].start;
				strongRoots[i].start = strongRoot->next;
				liftAll (strongRoot);
			}
		}
	}

	if (!strongRoots[0].start)
	{
		return NULL;
	}

	while (strongRoots[0].start)
	{
		strongRoot = strongRoots[0].start;
		strongRoots[0].start = strongRoot->next;
		strongRoot->label = 1;
		-- labelCount[0];
		++ labelCount[1];

		++ numRelabels;

		addToStrongBucket (strongRoot, &strongRoots[strongRoot->label]);
	}

	highestStrongLabel = 1;

	strongRoot = strongRoots[1].start;
	strongRoots[1].start = strongRoot->next;
	strongRoot->next = NULL;

	return strongRoot;
}



static void initializeRoot (Root *rt)
{
/*************************************************************************
initializeRoot
*************************************************************************/
	rt->start = NULL;
	rt->end = NULL;
}


static void initializeNode (Node *nd, const uint n)
{
/*************************************************************************
initializeNode
*************************************************************************/
	nd->label = 0;
	nd->excess = 0;
	nd->parent = NULL;
	nd->childList = NULL;
	nd->nextScan = NULL;
	nd->nextArc = 0;
	nd->numOutOfTree = 0;
	nd->arcToParent = NULL;
	nd->next = NULL;
	nd->visited = 0;
	nd->numAdjacent = 0;
	nd->number = n;
	nd->outOfTree = NULL;
}

static void initializeNodeSuper(NodeSuper *nd, const uint n)
{
	/*************************************************************************
	initializeNodeSuper
	*************************************************************************/
	nd->label = n;
	nd->numAdjacent = 0;
}

static void freeRoot (Root *rt)
{
/*************************************************************************
freeRoot
*************************************************************************/
	rt->start = NULL;
	rt->end = NULL;
}
static void freeMemory (void)
{
/*************************************************************************
freeMemory
*************************************************************************/
	uint i;

	for (i=0; i<numNodes; ++i)
	{
		freeRoot (&strongRoots[i]);
	}

	free(strongRoots);

	/*for (i=0; i<numNodes; ++i)
	{
		if (nodesList[i].outOfTree)
		{
			free(nodesList[i].outOfTree);
		}
	}*/

	free(nodesList);
	free(labelCount);
	free(arcList);
}
static void processRoot (Node *strongRoot)
{
/*************************************************************************
processRoot
*************************************************************************/
	Node *temp, *strongNode = strongRoot, *weakNode;
	Arc *out;

	strongRoot->nextScan = strongRoot->childList;

	if ((out = findWeakNode (strongRoot, &weakNode)))
	{
		merge (weakNode, strongNode, out);
		pushExcess (strongRoot);
		return;
	}

	checkChildren (strongRoot);

	while (strongNode)
	{
		while (strongNode->nextScan)
		{
			temp = strongNode->nextScan;
			strongNode->nextScan = strongNode->nextScan->next;
			strongNode = temp;
			strongNode->nextScan = strongNode->childList;

			if ((out = findWeakNode (strongNode, &weakNode)))
			{
				merge (weakNode, strongNode, out);
				pushExcess (strongRoot);
				return;
			}

			checkChildren (strongNode);
		}

		if ((strongNode = strongNode->parent))
		{
			checkChildren (strongNode);
		}
	}

	addToStrongBucket (strongRoot, &strongRoots[strongRoot->label]);
	++ highestStrongLabel;
}


static __inline void quickSort (Arc **arr, const uint first, const uint last)
{
/*************************************************************************
quickSort
*************************************************************************/
	int i=0, j, L=first, R=last, beg[MAX_LEVELS], end[MAX_LEVELS], temp=0;
	Arc *swap;

	if ((R-L) <= 5)
	{// Bubble sort if 5 elements or less
		for (i=R; (i>L); --i)
		{
			swap = NULL;
			for (j=L; j<i; ++j)
			{
				if (isExcess(arr[j]->flow - arr[j+1]->flow) < 0)//(arr[j]->flow < arr[j+1]->flow)
				{
					swap = arr[j];
					arr[j] = arr[j+1];
					arr[j+1] = swap;
				}
			}

			if (!swap)
			{
				return;
			}
		}

		return;
	}

	beg[0]=first;
	end[0]=last+1;

	while (i>=0)
	{
		L=beg[i];
		R=end[i]-1;

		if (L<R)
		{
			swap=arr[L];
			while (L<R)
			{
				while ((isExcess(arr[R]->flow - swap->flow)>=0) && (L<R)) //((arr[R]->flow >= swap->flow) && (L<R))
					R--;

				if (L<R)
				{
					arr[L]=arr[R];
					L++;
				}

				while ((isExcess(arr[L]->flow - swap->flow) <= 0) && (L<R)) //((arr[L]->flow <= swap->flow) && (L<R))
					L++;

				if (L<R)
				{
					arr[R]=arr[L];
					R--;
				}
			}

			arr[L] = swap;

			beg[i+1] = L+1;
			end[i+1] = end[i];
			end[i++] = L;

			if ((end[i]-beg[i]) > (end[i-1]-beg[i-1]))
			{
				temp=beg[i];
				beg[i]=beg[i-1];
				beg[i-1]=temp;
				temp=end[i];
				end[i]=end[i-1];
				end[i-1]=temp;
			}
		}
		else
		{
			i--;
		}
	}
}

static __inline void sort (Node * current)
{
/*************************************************************************
sort
*************************************************************************/
	if (current->numOutOfTree > 1)
	{
		quickSort (current->outOfTree, 0, (current->numOutOfTree-1));
	}
}

static __inline void minisort (Node *current)
{
/*************************************************************************
minisort
*************************************************************************/
	Arc *temp = current->outOfTree[current->nextArc];
	uint i, size = current->numOutOfTree;/*, tempflow = temp->flow;*/
	double tempflow = temp->flow;

	for(i=current->nextArc+1; ((i<size) && (isExcess(tempflow - current->outOfTree[i]->flow) < 0)); ++i)
	{
		current->outOfTree[i-1] = current->outOfTree[i];
	}
	current->outOfTree[i-1] = temp;
}


static /*ullint*/double checkOptimality (const uint gap)
{
/*************************************************************************
checkOptimality
*************************************************************************/
	uint i, check = 1;
	/*ullint*/double mincut = 0;
	/*llint*/ double *excess = NULL;

	Arc *tempArc;

	excess = (/*llint*/double *) malloc (numNodes * sizeof (double/*llint*/));
	if (!excess)
	{
		printf("Out of memory\n");
		exit(0);
	}

	// Pushing depicits from all sink adjacent nodes to the sink
	for (i=0; i<nodesList[sink-1].numOutOfTree; ++i)
	{
		tempArc = nodesList[sink-1].outOfTree[i];
		if (isExcess(tempArc->from->excess) < 0)
		{
			if (isExcess((tempArc->from->excess + /*(int)*/ tempArc->flow))  < 0)
			{
				// Excess is high enough to saturate the arc => Flow on residual arc is zeroed
				tempArc->from->excess += /*(int)*/ tempArc->flow;
				tempArc->flow = 0;
			}
			else
				// Excess is NOT high enough to saturate the arc => Excess is zeroed
			{
				tempArc->flow = /*(uint)*/ (tempArc->from->excess + /*(int)*/ tempArc->flow);
				tempArc->from->excess = 0;
			}
		}
	}

	for (i=0; i<numNodes; ++i)
	{
		excess[i] = 0;
	}

	for (i=0; i<numArcs; ++i)
	{
		if ((arcList[i].from->label >= gap) && (arcList[i].to->label < gap))
		{
			mincut += arcList[i].capacity;
		}

		if ((isExcess(arcList[i].flow - arcList[i].capacity)>0) || (isExcess(arcList[i].flow) < 0))
		{
			check = 0;
			printf("Warning - Capacity constraint violated on arc (%d, %d). Flow = %d, capacity = %d\n",
				   arcList[i].from->number,
				   arcList[i].to->number,
				   arcList[i].flow,
				   arcList[i].capacity);
		}
		excess[arcList[i].from->number - 1] -= arcList[i].flow;
		excess[arcList[i].to->number - 1] += arcList[i].flow;
	}

	for (i=0; i<numNodes; i++)
	{
		if ((i != (source-1)) && (i != (sink-1)))
		{
			if (isExcess(excess[i]))
			{
				check = 0;
				printf ("Warning - Flow balance constraint violated in node %d. Excess = %lld\n",
						i+1,
						excess[i]);
			}
		}
	}

	check = 1;

	if (isExcess(excess[sink-1] - mincut) != 0)//(excess[sink-1] != mincut)
	{
		check = 0;
		printf("Warning - Flow is not optimal - max flow does not equal min cut!\n");
	}

// 	if (check)
// 	{
// 		printf ("Solution checks as optimal. \t Max Flow: \t %lld\n", mincut);
// 	}

	free (excess);
	excess = NULL;
	return mincut;
}

static void decompose (Node *excessNode, const uint source, uint *iteration)
{
/*************************************************************************
decompose
*************************************************************************/
	Node *current = excessNode;
	Arc *tempArc;
	/*uint*/double bottleneck = excessNode->excess;

	// Find the bottleneck along a path to the source or on a cycle
	for ( ;(current->number != source) && (current->visited < (*iteration)) /*&& (current->nextArc < current->numOutOfTree)*/; // Added by Cheng
		 current = tempArc->from)
	{
		current->visited = (*iteration);
		tempArc = current->outOfTree[current->nextArc];

		if (isExcess(tempArc->flow - bottleneck) < 0) //(tempArc->flow < bottleneck)
		{
			bottleneck = tempArc->flow;
		}
	}

	if (current->number == source) // the DFS reached the source
	{
		excessNode->excess -= bottleneck;
		current = excessNode;

		// Push the excess all the way to the source
		while (current->number != source)
		{
			tempArc = current->outOfTree[current->nextArc]; // Pick arc going out of node to push excess to
			tempArc->flow -= bottleneck; // Push back bottleneck excess on this arc

			if (isFlow(tempArc->flow)) // If there is still flow on this arc do sort on this current node
			{
				minisort(current);
			} else {
				++ current->nextArc;
			}

			current = tempArc->from;
		}
		return;
	}

	++ (*iteration);

	// This part is temporarily added by Cheng Lu.
//	if (current->visited == (*iteration)-1)
		bottleneck = current->outOfTree[current->nextArc]->flow;
/*	else{
		current = excessNode;
		bottleneck = excessNode->excess;
	}*/

	while (current->visited < (*iteration)) //&& (current->nextArc < current->numOutOfTree)) // Added by Cheng
	{
		current->visited = (*iteration);
		tempArc = current->outOfTree[current->nextArc];

		if (isExcess(tempArc->flow - bottleneck) < 0)//(tempArc->flow < bottleneck)
		{
			bottleneck = tempArc->flow;
		}
		current = tempArc->from;
	}

	++ (*iteration);

	// This part is temporarily added by Cheng Lu.
/*	if (current->visited == (*iteration)-1)
		;
	else{
		current = excessNode;
	}*/

	while (current->visited < (*iteration)) //&& (current->nextArc < current->numOutOfTree)) // Added by Cheng
	{
		current->visited = (*iteration);

		tempArc = current->outOfTree[current->nextArc];
		tempArc->flow -= bottleneck;

		if (isFlow(tempArc->flow))
		{
			minisort(current);
			current = tempArc->from;
		} else {
			++ current->nextArc;
			current = tempArc->from;
		}
	}
}

static void recoverFlow (const uint gap)
{
/*************************************************************************
recoverFlow
*************************************************************************/
	uint iteration = 1;
	uint i, j;
	Arc *tempArc;
	Node *tempNode;

	Node **nodePtrArray;

	if ((nodePtrArray = (Node **) malloc (numNodes * sizeof (Node *))) == NULL)
	{
		printf("Out of memory\n");
		exit(0);
	}

	for (i=0; i < numNodes ; i++)
		nodePtrArray[i] = &nodesList[i];

	// Adding arcs FROM the Source to Source adjacent nodes.
	for (i=0; i<nodesList[source-1].numOutOfTree; ++i)
	{
		tempArc = nodesList[source-1].outOfTree[i];
		addOutOfTreeNode (tempArc->to, tempArc);
	}

	// Zeroing excess on source and sink nodes
	nodesList[source-1].excess = 0;
	nodesList[sink-1].excess = 0;

	for (i=0; i<numNodes; ++i)
	{
		tempNode = &nodesList[i];

		if ((i == (source-1)) || (i == (sink-1)))
		{
			continue;
		}

		if (tempNode->label >= gap) //tempNode is in SINK set
		{
			tempNode->nextArc = 0;
			if ((tempNode->parent) && (isFlow(tempNode->arcToParent->flow)))
			{
				addOutOfTreeNode (tempNode->arcToParent->to, tempNode->arcToParent);
			}

			for (j=0; j<tempNode->numOutOfTree; ++j)
			{ // go over all sink-set-node's arcs and look for arc with NO flow
				if (!isFlow(tempNode->outOfTree[j]->flow))
				{	// Remove arc with no flow
					-- tempNode->numOutOfTree;
					tempNode->outOfTree[j] = tempNode->outOfTree[tempNode->numOutOfTree];
					-- j;
				}
			}

			sort(tempNode);
		}
	}

	for (i=lowestPositiveExcessNode ; i < numNodes ; ++i)
	{
		tempNode = nodePtrArray[i];
		while (isExcess(tempNode->excess) > 0)
		{
			++ iteration;
			decompose(tempNode, source, &iteration);
		}
	}
	
	free(nodePtrArray);
}

static void readData(char *filename)
/*************************************************************************
readData
*************************************************************************/
{
	char buffer[32768];
	char lineIndicator;

	/* define parameters */
	uint numNodes;
	uint numArcs;
	uint arcCount = 0;
	uint currentNode;
	uint isSourceAssigned = 0;
	uint isSinkAssigned = 0;
	char sourceSinkIndicator;
	double lambdaLow;
	double lambdaHigh;
	uint from;
	uint to;
	double constantCapacity;
	double multiplierCapacity;

	// open input file
	FILE* f = fopen(filename, "r");
	assert(f);

	/* Read lines of input file */
	while (1)
	{
		if (fgets(buffer, sizeof buffer, f) != NULL)
		{
			switch (*buffer)
			{
			case 'p': /* initialize problem */
				sscanf(buffer, "p %u %u %lf %lf %u\n", &numNodes, &numArcs, &lambdaLow, &lambdaHigh, &roundNegativeWeight);
				printf("numNodes: %u, numArcs: %u, lambdaLow: %lf, lambdaHigh: %lf\n, round negative weight: %u\n", numNodes, numArcs, lambdaLow, lambdaHigh, roundNegativeWeight);

				if ((nodeListSuper = (NodeSuper *)malloc(numNodes * sizeof(NodeSuper))) == NULL)
				{
					printf("Could not allocate memory.\n");
					exit(0);
				}
				if ((arcListSuper = (ArcSuper *)malloc(numArcs * sizeof(ArcSuper))) == NULL)
				{
					printf("Could not allocate memory.\n");
					exit(0);
				}

				uint i;

				/* Initialization */
				for (i = 0; i < numNodes; ++i)
				{
					initializeNodeSuper(&nodeListSuper[i], i);
				}
				for (i = 0; i < numArcs; ++i)
				{
					initializeArcSuper(&arcListSuper[i]);
				}

				if (lambdaLow == lambdaHigh)
				{
					useParametricCut = 0;
				}

				break;
			case 'n':
				sscanf(buffer, "n %i %c\n", &currentNode, &sourceSinkIndicator);
				printf("Node: %u is of type %c\n", currentNode, sourceSinkIndicator);
				if (sourceSinkIndicator == 's')
				{
					/* check if source is valid */
					if (currentNode >= numNodes || currentNode < 0)
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
						source = currentNode;
						isSourceAssigned = 1;
						printf("Node %u is the source\n", source);
					}
				}
				else if (sourceSinkIndicator == 't')
				{
					/* check if sink is valid */
					if (currentNode >= numNodes || currentNode < 0)
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
						sink = currentNode;
						isSinkAssigned = 1;
						printf("Node %u is the sink\n", sink);
					}
				}
				else
				{
					printf("Node type: %c is unknown\n", sourceSinkIndicator);
					exit(0);
				}

				break;
			case 'a':
				sscanf(buffer, "a %u %u %lf %lf\n", &from, &to, &constantCapacity, &multiplierCapacity);
				printf("Arc: from: %u to %u with capacity: %lf and multiplier %lf\n", from, to, constantCapacity, multiplierCapacity);

				/* assign arc */
				if (from < 0 || to < 0 || from >= numNodes || to >= numNodes)
				{
					printf("Nodes are labeled from 0 to <number of nodes>  - 1\n");
					exit(0);
				}
				else if (from == to)
				{
					printf("Node %u has a self loop which is not allowed\n", from);
					exit(0);
				}
				else if (multiplierCapacity > 0 && from != source)
				{
					printf("Only source adjacent arcs can have a strictly positive capacity multiplier\n");
					exit(0);
				}
				else if (multiplierCapacity < 0 && to != sink)
				{
					printf("Only sink adjacent arcs can have a strictly negative capacity multiplier\n");
					exit(0);
				}
				else if (arcCount >= numArcs)
				{
					printf("Incorrect number of arcs specified\n");
					exit(0);
				}

				arcListSuper[arcCount].constant = constantCapacity;
				arcListSuper[arcCount].multiplier = multiplierCapacity;
				arcListSuper[arcCount].from = &nodeListSuper[from];
				arcListSuper[arcCount].to = &nodeListSuper[to];

				++arcCount;

				++nodeListSuper[from].numAdjacent;
				++nodeListSuper[to].numAdjacent;
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

	/* check if correct number of arcs has been specified */
	if (arcCount != numArcs)
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
}


//static void createMemStructures( )
//{
//	/* create memory structures */
//	for (i=0; i<numnodes; ++i)
//	{
//		createoutoftree (&nodeslist[i]);
//	}
//
//	for (i=0; i<numarcs; i++)
//	{
//		to = arclist[i].to->number;
//		from = arclist[i].from->number;
//		capacity = arclist[i].capacity;
//
//		if (!((source == to) || (sink == from) || (from == to)))
//		{
//			if ((source == from) && (to == sink))
//			{
//				arclist[i].flow = capacity;
//			} else if (to == sink) {
//				addoutoftreenode (&nodeslist[to-1], &arclist[i]);
//			} else {
//				addoutoftreenode (&nodeslist[from-1], &arclist[i]);
//			}
//		}
//	}
//}


static void pseudoflowPhase1 (void)
{
/*************************************************************************
pseudoflowPhase1
*************************************************************************/
	Node *strongRoot;
	while ((strongRoot = getHighestStrongRoot ()))
	{
		processRoot (strongRoot);
	}
}

static void printOutput (char *filename, double *times )
{
/*************************************************************************
printOutput
*************************************************************************/
	double stats[5];
	stats[0] = numArcScans;
	stats[1] = numMergers;
	stats[2] = numPushes;
	stats[3] = numRelabels;
	stats[4] = numGaps;
	
	// open outputFile
	FILE* f = fopen(filename,"w");

	// compute flow value
	double cut = 0;

	uint i;
	for (i = 0; i<numArcs; ++i)
	{
		if ((arcList[i].from->label >= numNodes) && (arcList[i].to->label < numNodes))
		{
			cut += arcList[i].capacity;
		}
	}

	// print value
	fprintf(f,"%lf\n", cut);

	// findSourceSet
	for (i=0; i<numNodes; ++i)
	{
		if (nodesList[i].label >= numNodes)
		{
			fprintf(f, "1");
		}
		else
		{
			fprintf(f, "0");
		}
		if (i < numNodes - 1)
			fprintf(f, ",");
		else
		{
			fprintf(f, "\n");
		}
	}

	for (i = 0; i < 3; i++)
	{
		fprintf(f, "%lf", times[i]);
		if (i < 2)
		{
			fprintf(f, ",");
		}
		else
		{
			fprintf(f, "\n");
		}
	}

	for (i = 0; i < 5; i++)
	{
		fprintf(f, "%lf", stats[i]);
		if (i < 4)
		{
			fprintf(f, ",");
		}
		else
		{
			fprintf(f, "\n");
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
	assert(argc == 3);

	/*ullint*/double flow = 0;
	double readStart, readEnd, initStart, initEnd, solveStart, solveEnd;

	numArcScans = 0;
	numMergers = 0;
	numPushes = 0;
	numRelabels = 0;
	numGaps = 0;

	readStart = clock();
	// readInput
	readData( argv[1] );
	readEnd = clock();

	exit(1);
	initStart = clock();
	// initialize
	simpleInitialization();
	initEnd = clock();

	solveStart = clock();
	// solve
	pseudoflowPhase1();
	solveEnd = clock();

	double times[3];
	times[0] = (readEnd - readStart)/CLOCKS_PER_SEC;
	times[1] = (initEnd - initStart)/CLOCKS_PER_SEC;
	times[2] = (solveEnd - solveStart)/CLOCKS_PER_SEC;

//	recoverFlow( numNodes );
//	flow = checkOptimality (numNodes);

	printOutput( argv[2], times);

	freeMemory ();

	return;

}