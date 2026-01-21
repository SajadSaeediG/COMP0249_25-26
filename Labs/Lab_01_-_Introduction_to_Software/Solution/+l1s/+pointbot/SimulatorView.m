classdef SimulatorView < ebe.graphics.EventGeneratorView

    properties(Access = protected)

        % The ground truth geometry
        groundTruthDrawer;

        % Obstacle geometry
        occluderGeometry;

        % Draw the observation mean and covariance
        gpsZRDrawer;

        % Geometry for the bearing measurements
        bearingMeasurementGeometry;

    end

    methods(Access = public)
        function obj = SimulatorView(config, eventGenerator)
            obj@ebe.graphics.EventGeneratorView(config, eventGenerator);
        end

        function start(obj)

            % Set up the drawing for the ground truth position
            xTrueColour = ebe.graphics.DistinguishableColours.assignColour('Ground Truth');
            obj.groundTruthDrawer = ebe.graphics.CircleDrawer(xTrueColour, 0.5);
            % Draw the map if it's defined
            scenario = obj.eventGenerator.scenario();

            if (isempty(scenario) == true)
                return
            end

            % Draw the different types of sensors if present
            if (isfield(scenario, 'sensors'))
                if (isfield(scenario.sensors, 'gps'))
                   obj.plotGPSAndOccluders(scenario.sensors.gps);
                end
                if (isfield(scenario.sensors, 'bearing'))
                    obj.plotBearingSensors(scenario.sensors.bearing);
                end
            end
        end

        function visualize(obj, events)

            % Update the ground truth position
            xTrue = obj.eventGenerator.xTrue();
            obj.groundTruthDrawer.update(xTrue([1;3]));

            % Now check all the events and process
            for e = 1 : length(events)                
                event = events{e};
                if (strcmp(event.type, 'gps') == true)
                    obj.visualizeGPSObservation(event);
                elseif (strcmp(event.type, 'bearing') == true)
                    obj.visualizeBearingObservation(event);
                end
            end
        end
    end

    methods(Access = protected)

        function plotGPSAndOccluders(obj, gps)

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
                    'EdgeColor', 'k', 'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'LineWidth', 2);
            end
        end

        function plotBearingSensors(obj, bearing)

            % Draw the bearing sensor observations
            zBearingColour = ebe.graphics.DistinguishableColours.assignColour('bearing');            
            obj.bearingMeasurementGeometry = plot(NaN, NaN, 'Color', zBearingColour, 'LineWidth', 2);

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

        function visualizeGPSObservation(obj, event)
            obj.gpsZRDrawer.update(event.data, event.covariance);
        end

        function visualizeBearingObservation(obj, event)
            scenario = obj.eventGenerator.scenario();
            numEvents = numel(event.info);
            bearing = event.data;
            xy = NaN(2, 3 * numEvents);
            for s = 1 : numel(event.info)
                sensor = scenario.sensors.bearing.sensors(event.info(s));
                xy(:,3*s-2) = sensor.position;
                theta = bearing(s) + deg2rad(sensor.orientation);
                xy(:,3*s-1) = sensor.position + sensor.detectionRange * [cos(theta); sin(theta)];
            end
            set(obj.bearingMeasurementGeometry, 'XData', xy(1, :), 'YData', xy(2, :))
        end
    end

end