classdef InitialPriorEdge < g2o.core.BaseUnaryEdge
    % InitialPriorEdge summary of InitialPriorEdge
    %
    % This class stores the factor representing the initial conditions
    % applied at time step 0.
    %
    % The measurement is directly of the full system state and has rank 3.

    % This edge stores the initial conditions
    
    methods(Access = public)
    
        function obj = InitialPriorEdge()
            % InitialPriorEdge for InitialPriorEdge
            %
            % Syntax:
            %   obj = VehicleKinematicsEdge()
            %
            % Description:
            %   Creates an instance of the InitialPriorEdge object.
            %   This predicts the state from one timestep to the next.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a InitialPriorEdge

            obj = obj@g2o.core.BaseUnaryEdge(3);
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

            obj.errorZ = obj.edgeVertices{1}.x - obj.z;
            obj.errorZ(3) = g2o.stuff.normalize_theta(obj.errorZ(3));
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

            obj.J{1} = eye(3);
        end        
    end
end