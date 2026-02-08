classdef VehicleKinematicsEdge < g2o.core.BaseBinaryEdge
    % VehicleKinematicsEdge summary of VehicleKinematicsEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is as follows.
    %
    % Define the rotation vector
    %
    %   M = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
    %
    % The new state is predicted from 
    %
    %   x_(k+1) = x_(k) + M * [vx;vy;theta]
    %
    % Note in this case the measurement is actually the mean of the process
    % noise. It has a value of 0. The error vector is given by
    %
    % e(x,z) = inv(M) * (x_(k+1) - x_(k))
    %
    % Note this requires estimates from two vertices - x_(k) and x_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k) and slot 2 contains x_(k+1).
    
    methods(Access = public)
        function obj = VehicleKinematicsEdge()
            % VehicleKinematicsEdge for VehicleKinematicsEdge
            %
            % Syntax:
            %   obj = VehicleKinematicsEdge()
            %
            % Description:
            %   Creates an instance of the VehicleKinematicsEdge object.
            %   This predicts the state from one timestep to the next.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a VehicleKinematicsEdge

            obj = obj@g2o.core.BaseBinaryEdge(3);            
        end

        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex. Note the error enters in a nonlinear manner, so the
            %   equation has to be rearranged to make the error the subject
            %   of the formulat
           
            % Rotation matrix from prior state
            priorX = obj.edgeVertices{1}.x;

            c = cos(priorX(3));
            s = sin(priorX(3));
            
            Mi = [c s 0;
                -s c 0;
                0 0 1];
            
            % Compute the error.
            obj.errorZ = Mi * (obj.edgeVertices{2}.x ...
                - priorX) - obj.z;
            
            % Wrap the heading error to -pi to pi
            obj.errorZ(3) = g2o.stuff.normalize_theta(obj.errorZ(3));
        end
        
        % Compute the Jacobians
        function linearizeOplus(obj)
            % LINEARIZEOPLUS Compute the Jacobians for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the Jacobians for the edge. Since we have two
            %   vertices which contribute to the edge, the Jacobians with
            %   respect to both of them must be computed.
            %

            priorX = obj.edgeVertices{1}.x;
            c = cos(priorX(3));
            s = sin(priorX(3));
            dx = obj.edgeVertices{2}.x - priorX;
            Mi = [c s 0;
                -s c 0;
                0 0 1];
            obj.J{2} = Mi;
            obj.J{1}(1, 1) = - c;
            obj.J{1}(1, 2) = - s;
            obj.J{1}(1, 3) = -dx(1) * s + dx(2) * c;
            obj.J{1}(2, 1) = s;
            obj.J{1}(2, 2) = - c;
            obj.J{1}(2, 3) = -dx(1) * c - dx(2) * s;
            obj.J{1}(3, 3) = -1;
        end

        function initialEstimate(obj)
            % INITIALESTIMATE Set the initial estimates of vertices.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial value of the predicted state using the
            %   prior state.

            priorX = obj.edgeVertices{1}.estimate();
            M = [cos(priorX(3)) -sin(priorX(3)) 0;
                sin(priorX(3)) cos(priorX(3)) 0;
                0 0 1];
            predictedX = priorX;
            predictedX = predictedX + M * obj.z;
            predictedX(3) = g2o.stuff.normalize_theta(predictedX(3));
            obj.edgeVertices{2}.setEstimate(predictedX);
        end

    end    
end