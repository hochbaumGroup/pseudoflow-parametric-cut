function [ lambdas, cuts ] = hpfCompleteParametric( lowProblem, highProblem )
% finds all breakpoints in lambda range ( low, high]. Note that if lambda =
% low is a breakpoint, then this is not included.

% tolerance for having equal slope
TOL = 1E-12;

display( lowProblem.lambdaValue )
display( highProblem.lambdaValue )

lowerSolved = false;
% solve lower bound problem if necessary
if isempty( lowProblem.optimalCutValue )
    lowerSolved = true;
    lowProblem.solve( true );
end

upperSolved = false;
% solve upper bound problem if nessecary
if isempty( highProblem.optimalCutValue )
    upperSolved = true;
    highProblem.solve( false );
end

lambdas = [];
cuts = zeros( lowProblem.nNodes, 0 );

% find intersection
if abs( highProblem.optimalCutLambda - lowProblem.optimalCutLambda ) > TOL
    lambdaIntersect =  ( lowProblem.optimalCutWeight - highProblem.optimalCutWeight ) / ( highProblem.optimalCutLambda - lowProblem.optimalCutLambda );
else
    lambdaIntersect = NaN;
end

% fprintf('Low: %f, High: %f, Intersect: %f\n', lowProblem.lambdaValue,highProblem.lambdaValue,lambdaIntersect )

if lambdaIntersect + TOL < highProblem.lambdaValue && lambdaIntersect - TOL > lowProblem.lambdaValue
    upperBoundIntersect = highProblem.copyNewLambda( lambdaIntersect );
    upperBoundIntersect.reduceUpper( highProblem );

    [ lambdasUpper, cutsUpper ] = hpfCompleteParametric( lowProblem, upperBoundIntersect );

    lowerBoundIntersect = lowProblem.copyNewLambda( lambdaIntersect );
    lowerBoundIntersect.reduceLower( lowProblem );

    [ lambdasLower, cutsLower ] = hpfCompleteParametric( lowerBoundIntersect, highProblem );

    % add new breakpoints
    lambdas = [ lambdas lambdasUpper lambdasLower ];
    cuts = [ cuts cutsUpper cutsLower ];
elseif lambdaIntersect == highProblem.lambdaValue
%     fprintf('Saving %f as breakpoint\n',highProblem.lambdaValue )
    lambdas = [ lambdas lambdaIntersect ];
    cuts = [ cuts lowProblem.optimalCut ];     
elseif lambdaIntersect == lowProblem.lambdaValue
%     fprintf('Saving %f as breakpoint\n',lowProblem.lambdaValue )
    lambdas = [ lambdas lambdaIntersect ];
    cuts = [ cuts lowProblem.optimalCut ];     
end

if lowerSolved && upperSolved
    cuts = [ cuts highProblem.optimalCut ];
    lambdas = [ lambdas highProblem.lambdaValue ];
end

[ lambdas, indices ] = unique( lambdas );
cuts = cuts( :, indices );