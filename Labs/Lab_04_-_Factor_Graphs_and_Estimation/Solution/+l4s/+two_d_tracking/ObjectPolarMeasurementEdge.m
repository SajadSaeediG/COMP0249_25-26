classdef ObjectPolarMeasurementEdge < g2o.core.BaseUnaryEdge
    % ObjectPolarMeasurementEdge summary of ObjectPolarMeasurementEdge
    %
    % This class stores an edge which represents the factor for observing
    % the vertex position using a range bearing sensors.
    %
    % The measurement model is
    %
    %    z_(k+1)=h[x_(k+1)]+w_(k+1)
    %
    % The measurements are r_(k+1) and beta_(k+1) and are given as follows.
    % The sensor is at (sx, sy, stheta).
    %
    %    dx = x(k+1) - sx; dy = y(k+1) - sy
    %
    %    r(k+1) = sqrt(dx^2+dy^2)
    %    beta(k+1) = atan2(dy, dx) - stheta
    %
    % The error term
    %    e(x,z) = z(k+1) - h[x(k+1)]
    %
    % However, remember that angle wrapping is required, so you will need
    % to handle this appropriately in compute error.

    properties(Access = protected)
        % The x,y and theta of the sensor
        sensorPose;
    end
    
    methods(Access = public)
    
        function obj = ObjectPolarMeasurementEdge()
            % ObjectPolarMeasurementEdge for ObjectPolarMeasurementEdge
            %
            % Syntax:
            %   obj = ObjectPolarMeasurementEdge()
            %
            % Description:
            %   Creates an instance of the ObjectPolarMeasurementEdge object.
            %   This is an observation of the particle's position.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ObjectPolarMeasurementEdge
            
            % Call the base constructor, specifying that this has a single
            % vertex and the dimension of the measurement is 2
            obj = obj@g2o.core.BaseUnaryEdge(2);

            % Set to a default value
            obj.sensorPose = zeros(3, 1);
        end
        
        function setSensorPose(obj, sensorPose)
            % SETSENSORPOSE Set the sensor pose.
            %
            % Syntax:
            %   obj.setSensorPose(sensorPose);
            %
            % Description:
            %   The range bearing sensor, sits in a given position and is
            %   oriented in a particular direction. This 
            %
            % Inputs:
            %   sensorPose - (3x1 double vector)
            %       The field are (x,y,theta) [theta in radians]

            obj.sensorPose = sensorPose;
        end
        
        function computeError(obj)
            % computeError Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the predicted and actual range-bearing measurement.
            %

            dXY = obj.edgeVertices{1}.x([1 3]) - obj.sensorPose(1:2);
            obj.errorZ(1) = obj.z(1) - norm(dXY);
            obj.errorZ(2) = g2o.stuff.normalize_theta(obj.z(2) - atan2(dXY(2), dXY(1)) + obj.sensorPose(3));
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
            
            dXY = obj.edgeVertices{1}.x([1 3]) - obj.sensorPose(1:2);
            r = norm(dXY);
            
            % The Jacobians are adapted from
            % https://github.com/petercorke/robotics-toolbox-matlab/blob/master/RangeBearingSensor.m 
            obj.J{1} = [-dXY(1)/r 0 -dXY(2)/r 0;
                dXY(2)/r^2 0 -dXY(1)/r^2 0];
        end        
    end
end