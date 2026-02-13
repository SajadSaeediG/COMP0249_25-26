classdef PlatformPredictionEdge < g2o.core.BaseBinaryEdge
    % PlatformPredictionEdge summary of PlatformPredictionEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is as follows.
    %
    % Define the rotation vector
    %
    %   M = dT * [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
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
    
    properties(Access = protected)
        % The length of the time step
        dT;
    end

    methods(Access = public)
        function obj = PlatformPredictionEdge(dT)
            % PlatformPredictionEdge for PlatformPredictionEdge
            %
            % Syntax:
            %   obj = PlatformPredictionEdge(dT);
            %
            % Description:
            %   Creates an instance of the PlatformPredictionEdge object.
            %   This predicts the state from one timestep to the next. The
            %   length of the prediction interval is dT.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a PlatformPredictionEdge

            assert(dT >= 0);
            obj = obj@g2o.core.BaseBinaryEdge(3);            
            obj.dT = dT;
        end
       
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of a platform.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the platform x_(k+1) given
            %   an estimate of the platform at time x_(k) and the control
            %   input u_(k+1)
            
            priorX = obj.edgeVertices{1}.estimate();
            u = obj.z;

            dt = obj.dT;
            
            % Rotation matrix scaled by timestep
            c = cos(priorX(3));
            s = sin(priorX(3));
            M = dt * [c -s 0;
                s  c 0;
                0  0 1];
            
            % Predict forward assuming no noise
            predictedX = priorX + M * u;
            predictedX(3) = g2o.stuff.normalize_theta(predictedX(3));
            obj.edgeVertices{2}.setEstimate(predictedX);
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
            
            priorX = obj.edgeVertices{1}.x;
            nextX = obj.edgeVertices{2}.x;
            
            c = cos(priorX(3));
            s = sin(priorX(3));

            dt = obj.dT;
            
            % Inverse of M (rotation scaled by dt)
            Mi = (1 / dt) * [c s 0;
                -s c 0;
                0 0 1];
            
            dx = nextX - priorX;
            
            % Compute the error
            obj.errorZ = Mi * dx - obj.z;
            
            % Wrap the heading error (in angle space) then convert to rate
            angleErr = g2o.stuff.normalize_theta(dx(3) - dt * obj.z(3));
            obj.errorZ(3) = angleErr / dt;
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
            nextX = obj.edgeVertices{2}.x;

            dt = obj.dT;
            
            c = cos(priorX(3));
            s = sin(priorX(3));
            dx = nextX - priorX;
            
            Mi = (1 / dt) * [c s 0;
                -s c 0;
                0 0 1];
            
            obj.J{2} = Mi;
            
            obj.J{1}(1, 1) = -c / dt;
            obj.J{1}(1, 2) = -s / dt;
            obj.J{1}(1, 3) = (-dx(1) * s + dx(2) * c) / dt;
            obj.J{1}(2, 1) = s / dt;
            obj.J{1}(2, 2) = -c / dt;
            obj.J{1}(2, 3) = (-dx(1) * c - dx(2) * s) / dt;
            obj.J{1}(3, 3) = -1 / dt;
        end
    end    
end
