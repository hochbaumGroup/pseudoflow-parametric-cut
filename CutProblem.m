classdef CutProblem < handle
    % CUTPROBLEM class for representing parametric minimum cut problems on
    % a directed graph. Problem assumes that the source-adjacent arcs are
    % non-decreasing in the parameter lambda whereas the sink-adjacent arcs
    % are non-increasing. Nodes are to be uniquely labeled in the set 1 :
    % nNodes, where nNodes is the number of nodes.
    % Author: Quico Spaen
    
    % TO DO:
    % - Add input checking to instance generation.
    
    properties
        % number of nodes in graph = n_1 + n_2 + n_3
        nNodes
        
        % 1 x n_1 row vector with labels of nodes in the source set 
        sourceSet
        % 1 x n_2 row vector with labels of nodes in the sink set
        sinkSet
        
        % 1 x n_3 row vector with labels of nodes that in neither source nor sink set.
        capLabels
        % n_3 x n_3 sparse matrix containing the capacities on the arc.
        % Element (i,j) corresponds with an arc from node capLabels[ i ] to
        % node capLabels[ j ].
        capacities
        
        % n_3 x 1 column vector containing the capacity of the source
        % adjacent arcs. The ith element is capacity of an arc from the
        % source to node capLabels[ i ].
        sourceWeights
        % n_3 x 1 column vector containing the capacity of the sink
        % adjacent arcs. The ith element is capacity of an arc from the
        % node capLabels[ i ] to the sink.
        sinkWeights
        
        % n_3 x 1 column vector with non-negative entries containing the 
        % lambda multiplier of the source adjacent arcs. 
        lambdaMultiplierSource
        % n_3 x 1 column vector with non-negative entries containing the 
        % lambda multiplier of the sink adjacent arcs. 
        lambdaMultiplierSink
        
        % scalar: weight of arc between source and sink
        sourceSinkWeight
        % scalar: lambda multiplier for arc from source to sink. Free in
        % sign.
        sourceSinkLambdaMultiplier
        
        % scalar: lambda value associated with this instance
        lambdaValue
        
        % nNodes x 1 vector. 1 indicates that a node is in the source set.
        % 0 indicates that the node is in the sink set.
        optimalCut
        % scalar: value of the optimal cut evaluated at lambdaValue
        optimalCutValue
        % scalar: weight of optimal cut function
        optimalCutWeight
        % scalar: lambda coefficient of optimal cut function
        optimalCutLambda
    end
    
    methods
        function obj = CutProblem( nNodes, sourceSet, sinkSet, capLabels, capacities, sourceWeights, sinkWeights, lambdaMultiplierSource, lambdaMultiplierSink, sourceSinkWeight, sourceSinkLambdaMultiplier, lambdaValue )
            % initialize cut instance.
            obj.nNodes = nNodes;
            obj.sourceSet = sourceSet;
            obj.sinkSet = sinkSet;
            obj.capLabels = capLabels;
            obj.capacities = capacities;
            obj.sourceWeights = sourceWeights;
            obj.sinkWeights = sinkWeights;
            obj.lambdaMultiplierSource = lambdaMultiplierSource;
            obj.lambdaMultiplierSink = lambdaMultiplierSink;
            obj.sourceSinkWeight = sourceSinkWeight;
            obj.sourceSinkLambdaMultiplier = sourceSinkLambdaMultiplier;
            obj.lambdaValue = lambdaValue;
        end
        
        function solve( obj, minimalTrue )
            lambdaSourceWeights = obj.sourceWeights + obj.lambdaValue .* obj.lambdaMultiplierSource;
            lambdaSinkWeights = obj.sinkWeights - obj.lambdaValue .* obj.lambdaMultiplierSink;
            lambdaSourceSinkWeight = obj.sourceSinkWeight + obj.lambdaValue * obj.sourceSinkLambdaMultiplier;
            
            nRemaining = length( obj.capLabels );
            capacitiesComplete = [ [ obj.capacities; lambdaSourceWeights'; zeros( 1, nRemaining ) ] ...
                zeros( nRemaining +2, 1 ), [ lambdaSinkWeights; lambdaSourceSinkWeight; 0 ] ];
            
            % inputs: Adjacency matrix graph, source, sink
            if minimalTrue
                [~,cut,~,~] = hpf( capacitiesComplete, nRemaining + 1, nRemaining + 2 );
            else
                [~,reversedCut,~,~] = hpf( capacitiesComplete', nRemaining + 2, nRemaining + 1 );
                cut = 1 - reversedCut;
            end
            
            remainingCut = logical( cut( 1 : nRemaining ) );
            obj.optimalCutWeight = sum( sum( obj.capacities( remainingCut, logical( 1 - remainingCut ) ) ) ) + sum( obj.sourceWeights( logical( 1 - remainingCut ) ) ) + sum( obj.sinkWeights( remainingCut ) ) + obj.sourceSinkWeight;
            obj.optimalCutLambda = sum( obj.lambdaMultiplierSource( logical( 1 - remainingCut ) ) ) - sum( obj.lambdaMultiplierSink( remainingCut ) ) + obj.sourceSinkLambdaMultiplier;
            
            % reconstruct original cut
            cutAll = zeros( obj.nNodes, 1);
            cutAll( obj.capLabels( remainingCut ) ) = 1;
            cutAll( obj.sourceSet ) = 1;
            
            obj.optimalCut = cutAll;
            obj.optimalCutValue = obj.optimalCutWeight + obj.lambdaValue * obj.optimalCutLambda;
        end
        
        function obj = copyNewLambda( self, lambdaValue )
            obj = CutProblem( self.nNodes, self.sourceSet, self.sinkSet, self.capLabels, self.capacities, self.sourceWeights, self.sinkWeights, self.lambdaMultiplierSource, self.lambdaMultiplierSink, self.sourceSinkWeight, self.sourceSinkLambdaMultiplier, lambdaValue );
        end
        
        function reduceUpper( self, highProblem )
            remainingIndicators = logical( highProblem.optimalCut( highProblem.capLabels ) );
            newSinkSetIndicators = not( remainingIndicators );
            
            self.sinkSet = [ self.sinkSet self.capLabels( newSinkSetIndicators ) ];
            self.capLabels = self.capLabels( remainingIndicators );
            
            self.sourceSinkWeight = self.sourceSinkWeight + sum( self.sourceWeights( newSinkSetIndicators ) );
            self.sourceSinkLambdaMultiplier = self.sourceSinkLambdaMultiplier + sum( self.lambdaMultiplierSource( newSinkSetIndicators ) );
            try
                self.sinkWeights = self.sinkWeights( remainingIndicators ) + sum( self.capacities( remainingIndicators, newSinkSetIndicators), 2 );
            catch
                self.sinkWeights = [];
            end
            self.sourceWeights = self.sourceWeights( remainingIndicators );
            self.capacities = self.capacities( remainingIndicators, remainingIndicators );
            self.lambdaMultiplierSource = self.lambdaMultiplierSource( remainingIndicators );
            self.lambdaMultiplierSink = self.lambdaMultiplierSink( remainingIndicators );
        end
        
        function reduceLower( self, lowProblem )
            newSourceSetIndicators = logical( lowProblem.optimalCut( lowProblem.capLabels ) );
            remainingIndicators = not( newSourceSetIndicators );
            
            self.sourceSet = [ self.sourceSet self.capLabels( newSourceSetIndicators ) ];
            self.capLabels = self.capLabels( remainingIndicators );
            
            self.sourceSinkWeight = self.sourceSinkWeight + sum( self.sinkWeights( newSourceSetIndicators ) );
            self.sourceSinkLambdaMultiplier = self.sourceSinkLambdaMultiplier - sum( self.lambdaMultiplierSink( newSourceSetIndicators ) );
            self.sinkWeights = self.sinkWeights( remainingIndicators );
            self.sourceWeights = self.sourceWeights( remainingIndicators )  + sum( self.capacities( newSourceSetIndicators, remainingIndicators ), 1 )';
            self.capacities = self.capacities( remainingIndicators, remainingIndicators );
            self.lambdaMultiplierSource = self.lambdaMultiplierSource( remainingIndicators );
            self.lambdaMultiplierSink = self.lambdaMultiplierSink( remainingIndicators );
        end
        
    end
    
end
