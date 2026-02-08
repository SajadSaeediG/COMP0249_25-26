classdef GPSMeasurementEdge < g2o.core.BaseUnaryEdge
    % GPSMeasurementEdge summary of GPSMeasurementEdge
    %
    % This class stores an edge which represents a GPS measurement. One
    % thing you can do is flip the signs of everything and everything will
    % still work. We do this here for illustration.
    %
    % The measurement is 2D and has the form
    %   z(1) - zx
    %   z(2) - zy
    %
    % We define the error vector as
    %
    %    e(x,z) = [x - zx; y - yz]
    
    methods(Access = public)
    
        function obj = GPSMeasurementEdge()
            % GPSMeasurementEdge for GPSMeasurementEdge
            %
            % Syntax:
            %   obj = GPSMeasurementEdge()
            %
            % Description:
            %   Creates an instance of the GPSMeasurementEdge object.
            %   This is an observation of the particle's position.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a GPSMeasurementEdge

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

            x = obj.edgeVertices{1}.estimate();
            obj.errorZ = obj.z - x(1:2);
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

            obj.J{1} = -eye(2, 3);
        end
    end
end