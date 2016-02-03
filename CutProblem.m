classdef CutProblem < handle
    
    properties
        nNodes
        
        sourceSet
        sinkSet
        
        capLabels
        capacities
        
        sourceWeights
        sinkWeights
        
        lambdaMultiplierSource
        lambdaMultiplierSink
        
        sourceSinkWeight
        sourceSinkLambdaMultiplier
        
        lambdaValue
        
        optimalCut
        optimalCutValue
        optimalCutWeight
        optimalCutLambda
        
    end
    
    methods
        function obj = CutProblem( nNodes, sourceSet, sinkSet, capLabels, capacities, sourceWeights, sinkWeights, lambdaMultiplierSource, lambdaMultiplierSink, sourceSinkWeight, sourceSinkLambdaMultiplier, lambdaValue )
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
