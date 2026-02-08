classdef PlatformStateVertex < g2o.core.BaseVertex
    % PlatformStateVertex summary of PlatformStateVertex
    %
    % This class stores the state of the vehicle. From the lab
    % description, the state has 3 dimensions, and the components are
    %
    % x(1) - x
    % x(2) - y
    % x(3) - theta
   
    properties(Access = protected)
        T;
    end
    
    methods(Access = public)
        function obj = PlatformStateVertex(T)
            % PlatformStateVertex for PlatformStateVertex
            %
            % Syntax:
            %   obj = VehicleStateVertex(time)
            %
            % Description:
            %   Creates an instance of the PlatformStateVertex object to
            %   store the estimate at time t.
            %
            % Inputs:
            %   time - (double)
            %       The time the vertex stores the estimate for.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a VehicleStateVertex
            
            obj=obj@g2o.core.BaseVertex(3);
            obj.T = T;
        end
        
        function T = time(obj)
            % TIME Return the time associated with this vertex
            %
            % Syntax:
            %   T = obj.time();
            %
            % Description:
            %   Return the time that this vertex is associated with.
            %
            %
            % Outputs:
            %   time - (double)
            %       The time for the vertex.

            T = obj.T;
        end
        
        function oplus(obj, update)
            % OPLUS Apply an incremental update to the state estimate.
            %
            % Syntax:
            %   obj.oplus(update);
            %
            % Description:
            %   The incremental update to the platform state estimate.
            %   Because this is (x,y,theta), we have to normalize the
            %   heading afterwards.
            %
            % Inputs:
            %   update - (3x1 double)
            %       Small perturbed update to the state estimate. This will
            %       have the same state dimension as the vertex state.
           
            % Add the update
            obj.x = obj.x + update;
            
            % Wrap the angle to [-pi,pi]
            obj.x(3) = g2o.stuff.normalize_theta(obj.x(3));
        end
    end
end