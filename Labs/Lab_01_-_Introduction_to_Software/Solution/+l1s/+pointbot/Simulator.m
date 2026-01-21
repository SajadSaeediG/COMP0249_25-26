classdef Simulator < ebe.core.EventBasedSimulator

    properties(Access = public)
        % The state
        x;

        % The map of sensors and occluders
        scenario;

        % The linear system used to predict the state and the observation
        systemModel;

        % The store of the ground truth state and observations
        xTrueStore;
        timeStore;

    end

    methods(Access = public)
        
        % Construct the object
        function obj = Simulator(config)
            obj@ebe.core.EventBasedSimulator(config);
            obj.systemModel = l1s.pointbot.SystemModel(config, config.perturbWithNoise);
        end
        
        % Get the ground truth
        function x = xTrue(obj)
            x = obj.x;
        end

        function [timeHistory, xTrueHistory] = history(obj)
            timeHistory = obj.timeStore;
            xTrueHistory = obj.xTrueStore;
        end

        function start(obj)

            % Call the base class start functions
            start@ebe.core.EventBasedSimulator(obj);

            % If a map is defined in the configuration file, load it
            if (isstruct(obj.config.scenario))
                obj.scenario = obj.config.scenario;
            end

            % Prestore the results
            obj.timeStore = NaN(1, obj.config.maximumStepNumber + 1);
            obj.xTrueStore = NaN(length(obj.config.x0), obj.config.maximumStepNumber + 1);

            % Set the initialization callbabk
            obj.eventGeneratorQueue.insert(0, @obj.initialize);

        end

        function scenario = getScenario(obj)
            scenario = obj.scenario;
        end

    end

    methods(Access = protected)

        function handlePredictForwards(obj, dT, time)
            obj.x = obj.systemModel.predictState(obj.x, dT);
        end

        function initialize(obj)
            assert(obj.stepNumber == 0)

            % Initialize the ground truth state
            obj.x = obj.config.x0 + obj.noiseScale * sqrtm(obj.config.P0) * randn(size(obj.config.x0));
            obj.initialized = true;

            % Construct and post the init event
            event = ebe.core.Event(obj.currentTime, 'init', obj.config.x0, obj.config.P0);
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);

            % Schedule the GPS measurement
            if (isfield(obj.config.scenario.sensors, 'gps') == true) && (obj.config.scenario.sensors.gps.enabled == true)
                obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.gps.measurementPeriod, ...
                    @obj.predictGPSObservation);
            end

            % Schedule the bearing measurement
            if (isfield(obj.config.scenario.sensors, 'bearing')  == true)  && (obj.config.scenario.sensors.bearing.enabled == true)
                obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.bearing.measurementPeriod, ...
                    @obj.predictBearingObservations);
            end

            % Schedule the timeout event
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.heartbeatPeriod, ...
                    @obj.generateHeartbeat);
        end

        function predictGPSObservation(obj)

            % We get the GPS observation only if the dotbot is not under an
            % occluder
            if (obj.isWithinOccluder() == true)
                event = ebe.core.Event(obj.currentTime, 'null_obs');
            else
                % Generate the observation
                [z, ~, R] = obj.systemModel.predictGPSObservation(obj.x);

                % Post the event
                event = ebe.core.Event(obj.currentTime, 'gps', z, R);
            end

            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.gps.measurementPeriod, @obj.predictGPSObservation);
        end

        function predictBearingObservations(obj)

            sensorIDs = obj.isDetectedByBearingSensors();

            if (isempty(sensorIDs) == true)
                event = ebe.core.Event(obj.currentTime, 'null_obs');
            else
                nz = numel(sensorIDs);
                z = zeros(1, nz);
                for s = 1 : nz
                    sensor = obj.scenario.sensors.bearing.sensors(sensorIDs(s));
                    z(s) = obj.systemModel.predictBearingObservation(obj.x, sensor.position, sensor.orientation);
                end
                event = ebe.core.Event(obj.currentTime, 'bearing', z, obj.config.scenario.sensors.bearing.sigmaR, sensorIDs); 
                if (nz > 1)
                    %keyboard
                end
            end

            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.bearing.measurementPeriod, @obj.predictBearingObservations);            
        end
        
        % Method to check if a point is within any occluder
        function isInside = isWithinOccluder(obj)

            if (isfield(obj.scenario.sensors.gps, 'occluders') == false)
                isInside = false;
                return
            end

            occluders = obj.scenario.sensors.gps.occluders;

            isInside = false;
            for i = 1:length(occluders)
                occluder = occluders(i);
                if obj.x(1) >= occluder.x_min && obj.x(1) <= occluder.x_max && ...
                   obj.x(3) >= occluder.y_min && obj.x(3) <= occluder.y_max
                    isInside = true;
                    return;
                end
            end
        end
        
        % Method to check if a point is visible to any sensor
        function [sensorIDs] = isDetectedByBearingSensors(obj)
            sensorIDs = [];
            sensors = obj.scenario.sensors.bearing.sensors;

            for s = 1:length(sensors)
                sensor = sensors(s);
                sensorPos = sensor.position;
                dx = obj.x(1) - sensorPos(1);
                dy = obj.x(3) - sensorPos(2);
                distance = hypot(dx, dy);
                
                if distance <= sensor.detectionRange
                    angleToPoint = atan2(dy, dx);
                    relativeAngle = angleToPoint - deg2rad(sensor.orientation);
                    relativeAngle = atan2(sin(relativeAngle), cos(relativeAngle));
                    
                    if abs(relativeAngle) <= deg2rad(sensor.detectionAngle)/2
                        sensorIDs(end + 1) = s;
                    end
                end
            end
        end

        function generateHeartbeat(obj)
            event = ebe.core.Event(obj.currentTime, 'heartbeat');
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.heartbeatPeriod, ...
                    @obj.generateHeartbeat);
        end

        function storeStepResults(obj)
            % Store
            obj.stepNumber;
            obj.timeStore(obj.stepNumber + 1) = obj.currentTime;
            obj.xTrueStore(:, obj.stepNumber + 1) = obj.x;
        end
    end
end