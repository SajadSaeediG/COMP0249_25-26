classdef CompassMeasurementEdge < g2o.core.BaseUnaryEdge
    % CompassMeasurementEdge summary of CompassMeasurementEdge
    %
    % This class stores an edge which represents a compass measurement.
    %
    % The measurement model is
    %
    %    z_(k+1) = theta_(k+1) + w_(k+1)
    %
    %
    % The error term
    %    e(x,z) = z(k+1) - theta_(k+1)
    %
    % However, remember that angle wrapping is required, so you will need
    % to handle this appropriately in compute error.

    methods(Access = public)
    
        function obj = CompassMeasurementEdge()
            % CompassMeasurementEdge for CompassMeasurementEdge
            %
            % Syntax:
            %   obj = CompassMeasurementEdge()
            %
            % Description:
            %   Creates an instance of the CompassMeasurementEdge object.
            %   This is an observation of the vehicle's heading angle
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a CompassMeasurementEdge.
            
            % Call the base constructor, specifying that this has a single
            % vertex and the dimension of the measurement is 2
            obj = obj@g2o.core.BaseUnaryEdge(1);
        end
        
        function computeError(obj)
            % computeError Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error.
            %
            x = obj.edgeVertices{1}.estimate();
            obj.errorZ = g2o.stuff.normalize_theta(x(3) - obj.z);
        end
        
        function linearizeOplus(obj)
            % linearizeOplus Compute the Jacobian of the error in the edge.
            %
            % Syntax:
            %   obj.linearizeOplus();
            %
            % Description:
            %   Compute the Jacobian of the error function with respect to
            %   the vertex.
            obj.J{1} = [0 0 1];
        end        
    end
end