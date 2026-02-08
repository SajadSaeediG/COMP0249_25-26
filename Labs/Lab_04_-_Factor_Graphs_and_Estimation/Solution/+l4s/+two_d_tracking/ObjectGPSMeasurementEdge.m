classdef ObjectGPSMeasurementEdge < g2o.core.BaseUnaryEdge
    % ObjectGPSMeasurementEdge summary of ObjectGPSMeasurementEdge
    %
    % This class stores the factor representing a direct measurement of the
    % object's position.
    % 
    % The measurement is 2D and has the form
    %   z(1) - zx
    %   z(2) - zy
    %
   
    methods(Access = public)
    
        function obj = ObjectGPSMeasurementEdge()
            % ObjectMeasurementEdge for ObjectGPSMeasurementEdge
            %
            % Syntax:
            %   obj = ObjectGPSMeasurementEdge()
            %
            % Description:
            %   Creates an instance of the ObjectGPSMeasurementEdge object.
            %   This is an observation of the particle's position.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ObjectGPSMeasurementEdge
            
            % Call the base constructor, specifying that this has a single
            % vertex and the dimension of the measurement is 2
            obj = obj@g2o.core.BaseUnaryEdge(2);
        end
        
        function computeError(obj)
            % computeError Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex.
            %
            % The error model e(x,z) has the form
            %   e(x, z) = [zx - x; zy -y]

            % Compute the error
            obj.errorZ = obj.z - obj.edgeVertices{1}.x([1 3]);
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
            %

            obj.J{1} = [-1 0 0 0;
                0 0 -1 0];
        end        
    end
end