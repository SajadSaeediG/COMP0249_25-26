% This edge describes the observation model
%
% z = x + w
%
% Specifically, the error is given by
%
% w = z -x
%
% Note the sign of the returned Jacobian.
%

classdef ParticleMeasurementEdge < g2o.core.BaseUnaryEdge
    % ParticleMeasurementEdge summary of ParticleMeasurementEdge
    %
    % This class stores an edge which represents the factor for observing
    % the vertex position. Specifically, suppose this edge has a
    % measurement m and is attached to a vertex v. The measurement model
    % has the form
    %
    % z = v.x + w
    %
    % where v.x is the value of x from the vertex, and w is the observation
    % noise. We rearrange this to get the error vector,
    %
    % w = e(z, v) = z - v.x
    %
   
    methods(Access = public)
    
        function obj = ParticleMeasurementEdge()
            % ParticleMeasurementEdge for ParticleMeasurementEdge
            %
            % Syntax:
            %   particleMeasurementEdge = ParticleMeasurementEdge()
            %
            % Description:
            %   Creates an instance of the ParticleMeasurementEdge object.
            %   This is an observation of the particle's x value.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ParticleMeasurementEdge
            
            % Call the base constructor, specifying that this has a single
            % vertex and the dimension of the measurement is 1
            obj = obj@g2o.core.BaseUnaryEdge(1);
        end
        
        function computeError(obj)
            % computeError Compute the error for the edge.
            %
            % Syntax:
            %   particleMeasurementEdge.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex.
            %
            %   In this case this has the form e(x, z) = z - v.x
            

            obj.errorZ = obj.z - obj.edgeVertices{1}.x;
        end
        
        function linearizeOplus(obj)
            % linearizeOplus Compute the Jacobian of the error in the edge.
            %
            % Syntax:
            %   particleMeasurementEdge.linearizeOplus();
            %
            % Description:
            %   Compute the Jacobian of the error function with respect to
            %   the vertex.
            %
            %   Since the error
            %
            %   In this case this has the form e(x, z) = z - v.x, the
            %   derivative with respect to v.x is -1.

            obj.J{1} = - 1;
        end        
    end
end