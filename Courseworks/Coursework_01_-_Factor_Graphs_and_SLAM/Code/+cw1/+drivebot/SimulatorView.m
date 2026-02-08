classdef SimulatorView < ebe.graphics.EventGeneratorView
    % SimulatorView summary of SimulatorView
    %
    % Visualizes the output of the trianglebot simulator

    properties(Access = protected)

        % The ground truth geometry
        groundTruthDrawer;

        % Obstacle geometry
        occluderGeometry;

        % Draw the observation mean and covariance
        gpsZRDrawer;

        % Geometry for the bearing measurements
        bearingObservationGeometry;

        slamObservationGeometry;

        % Information about the axes
        drawAxes;

        % If set to true, the plot is always centered on the ground truth
        % position
        centreAxes;

    end

    methods(Access = public)
        function obj = SimulatorView(config, eventGenerator)
            % SimulatorView Constructor for SimulatorView
            %
            % Syntax:
            %   slamSystemView = SimulatorView(config, eventGenerator)
            %
            % Description:
            %   Creates an instance of a SimulatorView object. Note that
            %   the simulator view will attach itself to whichever is the
            %   currently selected figure.
            %
            % Inputs:
            %   config - (struct)
            %       The configuration structure
            %   eventGenerator - (ebe.core.EventGenerator)
            %       The object to be visualized
            %
            % Outputs:
            %   slamSystemView - (handle)
            %       An instance of a SimulatorView

            obj@ebe.graphics.EventGeneratorView(config, eventGenerator);
            obj.drawAxes = axis;
            obj.centreAxes = false;
        end

        function setCentreAxesOnTruth(obj, centreAxes)
            % SETCENTREAXESONTRUTH Set flag to centre the axes
            %
            % Syntax:
            %   slamSystemView.setCentreAxesOnTruth(centreAxes)
            %
            % Description:
            %   This method sets an internal flag which governs whether the
            %   axes should be centred on the ground truth value of the
            %   platform or not.
            %
            % Inputs:
            %   centreAxes - (bool)
            %       Flag to centre the axes

            obj.centreAxes = centreAxes;
        end

        function start(obj)
            % START Start the simulator view
            %
            % Syntax:
            %   slamSystemView.start()
            %
            % Description:
            %   Draw static things which do not change (e.g., landmark
            %   maps, waypoints) and set up cached handles for graphics
            %   which are dynamically updated.

            % Allocate the ground truth colour first
            xTrueColour = ebe.graphics.DistinguishableColours.assignColour('Ground Truth');

            % Now allocate the colour for the waypoints
            waypointColour = ebe.graphics.DistinguishableColours.assignColour('waypoints');

            % Draw the map if it's defined
            scenario = obj.eventGenerator.scenario();

            % Draw the waypoints (these go on the bottom)
            if (isfield(obj.config, 'platformTrajectory'))
                waypoints = obj.config.platformTrajectory.waypoints;
                plot(waypoints(:, 1), waypoints(:, 2), 'o', 'Color', waypointColour, 'MarkerSize', 4, 'LineWidth', 2);

                % To draw the connected line, we have to add the first
                % point back again
                numWayPoints = size(waypoints, 1);
                plot(waypoints([1:numWayPoints 1], 1), waypoints([1:numWayPoints 1], 2), 'Color', waypointColour, 'LineWidth', 2);
            end

            if (isempty(scenario) == false)
                % Draw the landmarks if present
                if (isfield(scenario, 'landmarks'))
                   obj.plotLandmarks(obj.eventGenerator.landmarks());
                end
    
                % Draw the different types of sensors if present
                if (isfield(scenario, 'sensors'))
                    if (isfield(scenario.sensors, 'gps'))
                       obj.plotGPSAndOccluders(scenario.sensors.gps);
                    end
                    if (isfield(scenario.sensors, 'bearing'))
                        obj.plotBearingSensors(scenario.sensors.bearing);
                    end
                    if (isfield(scenario.sensors, 'slam'))
                        obj.plotSLAMSensor(scenario.sensors.slam);
                    end
                end                
            end



            % Set up the drawing for the ground truth position
            obj.groundTruthDrawer = ebe.graphics.OrientedTriangleDrawer(xTrueColour);

            % Now 
        end

        function [handle, entry] = legendEntries(obj)
            xTrueColour = ebe.graphics.DistinguishableColours.assignColour('Ground Truth');
            handle = plot(NaN, NaN, 'Color', xTrueColour, 'LineWidth', 2);
            entry = 'True';
        end


        function visualize(obj, events)
            % VISUALIZE Update the visualizations.
            %
            % Syntax:
            %   slamSystemView.visualize(events);
            %
            % Description:
            %   Updates the graphical representation. This includes showing
            %   data from the simulator such as ground truth position and
            %   recently generated events
            %
            % Inputs:
            %   events - (cell array of ebe.core.Events)
            %       The latest list of events to visualize

            % Update the ground truth position
            xTrue = obj.eventGenerator.xTrue();
            obj.groundTruthDrawer.update(xTrue);

            % Now check all the events and process
            for e = 1 : length(events)                
                event = events{e};
                if (strcmp(event.type, 'gps') == true)
                    obj.visualizeGPSObservation(event);
                elseif (strcmp(event.type, 'bearing') == true)
                    obj.visualizeBearingObservation(event);
                elseif (strcmp(event.type, 'slam') == true)
                    obj.visualizeSLAMObservations(event);
                end
            end

            % Reposition the axes if necessary
            if (obj.centreAxes == true)
                axes = obj.drawAxes();
                axes(1:2) = axes(1:2) + xTrue(1);
                axes(3:4) = axes(3:4) + xTrue(2);
                axis(axes);
            end
        end
    end

    methods(Access = protected)

        function plotGPSAndOccluders(obj, gps)
            % PLOTGPSANDOCCLUDERS Set up graphics to plot the GPS sensor.
            %
            % Syntax:
            %   slamSystemView.plotGPSAndOccluders(gps)
            %
            % Description:
            %   Set up the graphics for the GPS, including the occluders
            %   and the measurement itself.
            %
            % Inputs:
            %   gps - (struct)
            %       The part of the system configuration which specifies
            %       the GPS sensor properties.

            % Draw the GPS observations
            zGPSColour = ebe.graphics.DistinguishableColours.assignColour('gps');            
            obj.gpsZRDrawer = ebe.graphics.MeanCovarianceDrawer(zGPSColour);

            if (isfield(gps, 'occluders') == false)
                return
            end

            numOccluders = numel(gps.occluders);

            obj.occluderGeometry = cell(numOccluders, 1);
            for i = 1:numOccluders
                occluder = gps.occluders(i);
                obj.occluderGeometry{i}=rectangle('Position', [occluder.x_min, occluder.y_min, ...
                    occluder.x_max - occluder.x_min, occluder.y_max - occluder.y_min], ...
                    'EdgeColor', 'k', 'FaceColor', [0.8, 0.8, 0.8, 0.7], 'LineWidth', 2);
            end
        end

        function plotBearingSensors(obj, bearing)
            % PLOTBEARINGSSENSORS Set up graphics to plot the bearing
            % sensor.
            %
            % Syntax:
            %   slamSystemView.plotBearingSensors(bearing)
            %
            % Description:
            %   Set up the graphics for the bearing sensors. This includes
            %   the detection cones.
            %
            % Inputs:
            %   bearing - (struct)
            %       The part of the system configuration which specifies
            %       the bearing sensor properties.

            zBearingColour = ebe.graphics.DistinguishableColours.assignColour('bearing');            
            obj.bearingObservationGeometry = plot(NaN, NaN, 'Color', zBearingColour);

            for s = 1 : numel(bearing.sensors)
                sensor = bearing.sensors(s);

                % Extract sensor properties
                pos = sensor.position;
                orientation = sensor.orientation;
                range = sensor.detectionRange;
                angle = sensor.detectionAngle;

                % Compute circular arc points
                theta = linspace(orientation - angle / 2, orientation + angle / 2, 100) * pi / 180;
                arcX = pos(1) + range * cos(theta);
                arcY = pos(2) + range * sin(theta);

                % Complete the wedge by connecting to the sensor's position
                wedgeX = [pos(1), arcX, pos(1)];
                wedgeY = [pos(2), arcY, pos(2)];

                % Plot sensor wedge with transparency
                fill(wedgeX, wedgeY, [1, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

                % Plot edges of the wedge
                line([pos(1), arcX(1)], [pos(2), arcY(1)], 'Color', 'r', 'LineWidth', 2);
                line([pos(1), arcX(end)], [pos(2), arcY(end)], 'Color', 'r', 'LineWidth', 2);

                % Plot circular arc
                plot(arcX, arcY, 'r', 'LineWidth', 2);
            end
        end

        function plotSLAMSensor(obj, ~)
            % PLOTSLAMSENSOR Set up graphics to plot the SLAM sensor
            % sensor.
            %
            % Syntax:
            %   slamSystemView.plotSLAMSensor()
            %
            % Description:
            %   Set up the graphics for the bearing sensors. This includes
            %   the detection cones.
            %
            % Inputs:
            %   bearing - (struct)
            %       The part of the system configuration which specifies
            %       the bearing sensor properties.

            slamSensorColour = ebe.graphics.DistinguishableColours.assignColour('slam_sensor_colour');
            obj.slamObservationGeometry = plot(NaN, NaN, 'Color', slamSensorColour, 'LineWidth', 2);
        end

        function plotLandmarks(obj, landmarks)
            % PLOTLANDMARKS Set up graphics to plot the landmarks
            % identified by SLAM
            %
            % Syntax:
            %   slamSystemView.plotLandmarks(landmarks)
            %
            % Description:
            %   Set up the graphics for the landmarks sensors.
            %
            % Inputs:
            %   landmarks - (struct)
            %       The part of the system configuration which specifies
            %       the landmark properties.

            landmarkColour = ebe.graphics.DistinguishableColours.assignColour('slam_true_landmarks');
            plot(landmarks(1, :), landmarks(2, :), '+', 'Color', landmarkColour, 'MarkerSize', 7, 'LineWidth', 2);
        end

        function visualizeGPSObservation(obj, event)
            % VISUALIZEGPSOBSERVATION Update the graphics for the GPS
            % observation
            %
            % Syntax:
            %   slamSystemView.visualizeGPSObservation(event)
            %
            % Description:
            %   Update the graphics (mean and covariance ellipse) for the
            %   GPS observation, based on the data stored in the event.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       A GPS observation event.

            obj.gpsZRDrawer.update(event.data, event.covariance);
        end

        function visualizeBearingObservation(obj, event)
            % VISUALIZEBEARINGOBSERVATION Update the graphics for the
            % bearing observation.
            %
            % Syntax:
            %   slamSystemView.visualizeBearingObservation(bearing)
            %
            % Description:
            %   Update the measurement (rays) in all the bearing sensors
            %   that observe the platform using data stored in the event.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       A bearing observation event.

            map = obj.eventGenerator.map();
            numObservations = numel(event.info);
            bearing = event.data;
            xy = NaN(2, 3 * numObservations);
            for s = 1 : numObservations
                sensor = map.sensors.bearing.sensors(event.info(s));
                xy(:,3*s-2) = sensor.position;
                theta = bearing(s) + deg2rad(sensor.orientation);
                xy(:,3*s-1) = sensor.position + sensor.detectionRange * [cos(theta); sin(theta)];
            end
            set(obj.bearingObservationGeometry, 'XData', xy(1, :), 'YData', xy(2, :))
        end

        function visualizeSLAMObservations(obj, event)
            % VISUALIZESLAMOBSERVATIONS Update the graphics for the
            % SLAM observations.
            %
            % Syntax:
            %   slamSystemView.visualizeSLAMObservations(event)
            %
            % Description:
            %   Draw the rays which show all the SLAM landmarks observed at
            %   a given time, using data from the event.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The landmark observation event.

            xTrue = obj.eventGenerator.xTrue();

            z = event.data;

            numObservations = numel(event.info);
            xy = NaN(2, 3 * numObservations);
            for s = 1 : numObservations
                beta = z(2,s) + xTrue(3);
                xy(:, 3*s-2) = xTrue(1:2);
                xy(:, 3*s-1) = xy(:, 3*s-2) + z(1, s) * [cos(beta);sin(beta)];
            end
            set(obj.slamObservationGeometry, 'XData', xy(1, :), 'YData', xy(2, :));
        end

    end

end