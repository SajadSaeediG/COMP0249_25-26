classdef ObjectProcessModelEdge < g2o.core.BaseBinaryEdge
    % ObjectProcessModelEdge summary of ObjectProcessModelEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is
    %    x_(k+1)=F_(k+1)*x_(k)+v_(k+1)
    %
    % Note in this case the measurement is actually the mean of the process
    % noise. It has a value of 0. The error vector
    %
    % e(x,z) = v_(k+1) = x_(k+1) - F_(k+1)*x_(k)
    %
    % Note this requires estimates from two vertices - x_(k) and x_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k) and slot 2 contains x_(k+1).

    properties(Access = protected)
        
        % The state transition matrix
        F;
        
    end
    
    methods(Access = public)
        function this = ObjectProcessModelEdge()
            % ObjectProcessModelEdge for ObjectProcessModelEdge
            %
            % Syntax:
            %   obj = ObjectProcessModelEdge()
            %
            % Description:
            %   Creates an instance of the ObjectProcessModelEdge object.
            %   This predicts the state from one timestep to the next.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ObjectProcessModelEdge
            
            this = this@g2o.core.BaseBinaryEdge(4);   
            
            % Create the state transtion matrix. Set it to NaN to trigger
            % errors if we don't assign a value to it.
            this.F = NaN(4, 4);            
        end
        
        function setF(this, F)
            % SETF Set the F matrix on this edge
            %
            % Syntax:
            %   obj.setF(F);
            %
            % Description:
            %   Set the matrix used in the linear model. We specify it here
            %   outside the edge so that the controlling programme can feed
            %   it in. Note that the process noise needs to match the form
            %   of F. However, the process noise is set using
            %   setInformation, and that is handled by another call.
            %
            % Inputs:
            %   F - (4x4 double matrix)
            %       The state transition matrix.

            this.F = F;
        end
       
        function computeError(this)
            % COMPUTEERROR Compute the error for the edge.
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
            %   e(x, z) = x_(k+1) - F * x_(k)

            this.errorZ = this.edgeVertices{2}.x - ...
                this.F * this.edgeVertices{1}.x;
        end
        
        function linearizeOplus(this)
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
            this.J{1} =  - this.F;
            this.J{2} = eye(4);
        end        
        
    end
    
end