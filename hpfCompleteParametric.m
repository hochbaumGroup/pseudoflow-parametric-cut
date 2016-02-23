function [ lambdas, cuts ] = hpfCompleteParametric( lowProblem, highProblem )
% HPFCOMPLETEPARAMETRIC solves a parametric s-t minimum cut problem by solving the maximum flow problem with Hochbaum's Pseudoflow.
% The algorithm finds all breakpoints for which the cut changes as a function of lambda in the range
% [ lambdaLow, lambdaHigh] by recursively concluding that the interval contains 0, 1, or more breakpoints.
% If the interval contains more than 1 breakpoint, then the interval is
% split into two interval, each of which contains a breakpoint.
%
% Parametric cut/flow problems allow for a linear function with input lambda
% on source or sink adjacent arcs. Arcs that are adjacent to source should
% be non-decreasing in lambda and sink adjacent arcs should be
% non-increasing in lambda. The algorithm is able to deal with the reverse
% configuration (non-increasing on source adjacent arcs and non-decreasing
% on sink adjacent arcs) by flipping source and sink and reversing the
% direction of the arcs.
%
% Theoretical background:
% - D. Hochbaum (2008). The Pseudoflow Algorithm: A New Algorithm for the 
%   Maximum-Flow Problem. Operations Research 56(4):992-1009.
% - D. Hochbaum, (2003) Efficient Algorithms for the Inverse Spanning-Tree 
%   Problem. Operations Research 51(5):785-797.
%
% ---------- INPUT ---------------
% lowProblem: minimum cut problem instance instantiated as an object of
%       class Cutproblem. Should have lambda value of lambdaLow.
% highProblem: same minimum cut problem instance but with lambda value
%       equal to lambdaHigh.
% ---------- OUTPUT --------------
% lambdas: 1 x k sorted vector containing all k-1 breakpoints and adds the
%       value lambdaHigh.
% cuts: n x k matrix containing the cuts corresponding to the lambdas
%       vector. The jth column of the cuts matrix is optimal for lambda in the
%       range ( lambdas[j-1], lambdas[j] ]. If the (i,j) element of the
%       cuts matrix is 1, then node i is in the source set of the cut for 
%       lambda in ( lambdas[j-1], lambdas[j] ]. If the element is 0, then
%       the node is in the sink set of the cut.
% ---------------------------------------------------------------------
% Author: Quico Spaen
% ---------------------------------------------------------------------


% tolerance for having equal slope
TOL = 1E-12;

display( lowProblem.lambdaValue )
display( highProblem.lambdaValue )


lowerSolved = false;
% solve lower bound problem if necessary by finding minimal source set
if isempty( lowProblem.optimalCutValue )
    lowerSolved = true;
    lowProblem.solve( true );
end

upperSolved = false;
% solve upper bound problem if nessecary by finding maximal source set
if isempty( highProblem.optimalCutValue )
    upperSolved = true;
    highProblem.solve( false );
end

% initialize output variables
lambdas = [];
cuts = zeros( lowProblem.nNodes, 0 );

% find lambda value for which the optimal cut functions (expressed as a 
% function of lambda) for the lower bound and upper bound problem intersect.
if abs( highProblem.optimalCutLambda - lowProblem.optimalCutLambda ) > TOL
    lambdaIntersect =  ( lowProblem.optimalCutWeight - highProblem.optimalCutWeight ) / ( highProblem.optimalCutLambda - lowProblem.optimalCutLambda );
else % conclude that there is no intersection if denominator is too close to zero.
    lambdaIntersect = NaN;
end

% check cases depending on intersection value of the lower bound and upper
% bound optimal cut.
if lambdaIntersect + TOL < highProblem.lambdaValue && lambdaIntersect - TOL > lowProblem.lambdaValue
    % if intersection occurs strictly within the interval, then there are
    % at least 2 breakpoints. lambdaIntersect is guaranteed to separate the
    % interval into subintervals each containing at least 1 breakpoint.
    
    % Create new instance of upper bound problem with contracted sink set 
    % and lambda value equal to lambda intersect. The nodes that are in the 
    % sinkset for lambdaHigh are guaranteed to be in the sink set for
    % lambda <= lambdaIntersect.
    upperBoundIntersect = highProblem.copyNewLambda( lambdaIntersect );
    upperBoundIntersect.reduceUpper( highProblem );

    % recurse for lower subinterval
    [ lambdasUpper, cutsUpper ] = hpfCompleteParametric( lowProblem, upperBoundIntersect );

    % Create new instance of lower bound problem with contracted source set 
    % and lambda value equal to lambda intersect. The nodes that are in the 
    % source set for lambdaLow are guaranteed to be in the source set for
    % lambda >= lambdaIntersect.
    lowerBoundIntersect = lowProblem.copyNewLambda( lambdaIntersect );
    lowerBoundIntersect.reduceLower( lowProblem );

    % recurse for higher subinterval
    [ lambdasLower, cutsLower ] = hpfCompleteParametric( lowerBoundIntersect, highProblem );

    % combine lambdas and cuts from subintervals
    lambdas = [ lambdas lambdasUpper lambdasLower ];
    cuts = [ cuts cutsUpper cutsLower ];
    
elseif lambdaIntersect == highProblem.lambdaValue
    % if lambda intersect is equal to upper bound, then lambdaHigh is a
    % breakpoint and no further recursion necessary.
    lambdas = [ lambdas lambdaIntersect ];
    cuts = [ cuts lowProblem.optimalCut ];
    
elseif lambdaIntersect == lowProblem.lambdaValue
	% if lambda intersect is equal to lower bound, then lambdaLow is a
    % breakpoint and no further recursion necessary.
    lambdas = [ lambdas lambdaIntersect ];
    cuts = [ cuts lowProblem.optimalCut ];     
end

% add cut corresponding to lambdaHigh iff first recursion level
if lowerSolved && upperSolved
    cuts = [ cuts highProblem.optimalCut ];
    lambdas = [ lambdas highProblem.lambdaValue ];
end

% lambdas will appear twice. Once as lambdaLow and once as lambdaHigh. This
% removes duplicates.
[ lambdas, indices ] = unique( lambdas );
cuts = cuts( :, indices );