classdef WaypointController < ebe.core.ConfigurableComponent
    % WaypointController summary of WaypointController
    %
    % The controller guides the movement of the trianglebot so it drives
    % close to a set of waypoints. The code here uses a driven-steered
    % vehicle model developed by Tim Bailey at ACFR. The control signals
    % are internally converted to the right form which can be used by
    % trianglebot.

    properties(Access = protected)

        % An array of (2x1) waypoints
        waypoints;

        % The number of waypoints
        numWaypoints;

        % Which waypoint is the robot heading towards now?
        waypointIndex;

        % The control input.         
        u;

        % Configuration of the platform
        platformConfig;

        % Configuration of the controller
        controllerConfig;

        % Flag to indicate if the robot should repeat the waypoint sequence
        repeatVisitingWaypoints;
        
        % Distance to waypoint
        maxAcceptableDistanceFromWaypoint;

    end

    methods(Access = public)

        function obj = WaypointController(config)
            % WaypointController Constructor for WaypointController
            %
            % Syntax:
            %   waypointController = WaypointController(config)
            %
            % Description:
            %   Creates an instance of a WaypointController object.
            %
            % Inputs:
            %   config - (struct)
            %       The configuration structure
            %
            % Outputs:
            %   waypointController - (handle)
            %       An instance of a WaypointController

            obj@ebe.core.ConfigurableComponent(config);
            obj.platformConfig = obj.config.platform;
            obj.controllerConfig = obj.config.platform.controller;
    
            if (isfield(obj.config.platformTrajectory, 'repeat') == true)
                obj.repeatVisitingWaypoints = obj.config.platformTrajectory.repeat;
            else
                obj.repeatVisitingWaypoints = true;
            end

            obj.maxAcceptableDistanceFromWaypoint = 1;
        end

        function start(obj)
            % START start the waypoint controller
            %
            % Syntax:
            %   waypointController.start()
            %
            % Description:
            %   Sets the controller up so that it's ready to go. This
            %   includes extracting the waypoints from the configuration
            %   file, and setting the first waypoint as the target.

            obj.waypoints = obj.config.platformTrajectory.waypoints';
            obj.numWaypoints = size(obj.waypoints, 2);
            obj.waypointIndex = min(1, obj.numWaypoints);
            obj.u=[0;0];

        end

        function u = computeControlInputs(obj, x)
            % COMPUTECONTROLINPUTS compute the next output from the
            % controller
            %
            % Syntax:
            %   u = waypointController.computeControlInputs(x);
            %
            % Description:
            %   Computes the control input which will steer the trianglebot
            %   towards the next waypoint. If the robot draws sufficiently
            %   close to the waypoint, the next waypoint is selected. If
            %   the waypoint list is completed, the platform can optionally
            %   go back to the first waypoint.
            %
            %   The controller itself is for a driven steered vehicle and
            %   has maximum threshold rates.
            %
            % Outputs:
            %   u - (2x1 double)
            %       The control input (speed, angular velocity) or an empty
            %       matrix if the controller has run out of waypoints and
            %       repeat is not set.
            
            % Work out distance to the target waypoint
            dX = obj.waypoints(:, obj.waypointIndex) - x(1:2);
            d = norm(dX);
            
            % If sufficiently close, we can look at switching to the next
            % waypoing
            if (d < obj.maxAcceptableDistanceFromWaypoint)

                % Handle the case if we have now visited all waypoints
                if (obj.waypointIndex == obj.numWaypoints)
                    if (obj.repeatVisitingWaypoints == true)
                        obj.waypointIndex = 0;
                    else
                        u = [];
                        return;
                    end
                end
    
                % Increment the new waypoint
                obj.waypointIndex = obj.waypointIndex + 1;
                
                % Switch to computing errors from the new waypoint
                dX = obj.waypoints(:, obj.waypointIndex) - x(1:2);
                d = norm(dX);
            end
            
            % Compute the speed. We first clamp the acceleration, and then
            % clamp the maximum and minimum speed values.
            diffSpeed = 0.1 * d - obj.u(1);
            maxDiffSpeed = obj.controllerConfig.maxAcceleration * obj.controllerConfig.odomUpdatePeriod;
            diffSpeed = min(maxDiffSpeed, max(-maxDiffSpeed, diffSpeed));
            obj.u(1) = max(obj.controllerConfig.minSpeed, min(obj.controllerConfig.maxSpeed, obj.u(1) + diffSpeed));

            % Compute the steer angle. We first clamp the rate of change,
            % and then clamp the maximum and minimum steer angles.
            diffDelta = g2o.stuff.normalize_theta(atan2(dX(2), dX(1)) - x(3) - obj.u(2));
            maxDiffDelta = obj.controllerConfig.maxDiffDeltaRate * obj.controllerConfig.odomUpdatePeriod;
            diffDelta = min(maxDiffDelta, max(-maxDiffDelta, diffDelta));
            obj.u(2) = min(obj.controllerConfig.maxDelta, max(-obj.controllerConfig.maxDelta, obj.u(2) + diffDelta));

            % Now work out what the control input vector. The main thing is
            % that we have to compute the platform angular velocity from
            % the bicycle model.
            psiDot = obj.u(1) * sin(obj.u(2)) / obj.platformConfig.B;
            u = [obj.u(1); 0; psiDot];
        end
    end
end