classdef VehicleStateVertex < g2o.core.BaseVertex
    % VehicleStateVertex summary of VehicleStateVertex
    %
    % This class stores the state of the vehicle. From the lab
    % description, the state has 3 dimensions, and the components are
    %
    % x(1) - x
    % x(2) - y
    % x(3) - theta
   
    methods(Access = public)
        function obj = VehicleStateVertex()
            % VehicleStateVertex for VehicleStateVertex
            %
            % Syntax:
            %   obj = VehicleStateVertex()
            %
            % Description:
            %   Creates an instance of the VehicleStateVertex object. This
            %   is four dimensional Euclidean state, and so no special
            %   oplus operator is required.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a VehicleStateVertex
            
            % Call the base class constructor, specifying the dimension of
            % the state.
            obj=obj@g2o.core.BaseVertex(3);
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