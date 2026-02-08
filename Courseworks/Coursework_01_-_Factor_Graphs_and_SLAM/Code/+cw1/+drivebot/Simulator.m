classdef Simulator < ebe.core.EventBasedSimulator
    % Simulator summary of Simulator
    %
    % Simulates the scenario for this exercise. The scenario consists of
    % the drivebot, following a set of waypoints, through an environment
    % populated by a set of landmarks. Depending upon the configuration
    % used, different trajectories, landmarks, sensing systems and noise
    % levels are available.
    %
    % In this particular example, the ground truth platform does not
    % experience process noise. Rather, control noisy values of control
    % inputs are observed.

    properties(Access = public)
        % The platform state
        x;

        % The latest control input
        u;

        % The covariance in the control input measurement
        sigmaU;

        % The standard deviation in the control input measurement, used to
        % simulate noise.
        sigmaUSqrtm;

        % The scenario, which describes things like landmarks and sensors
        scenario;

        % The landmarks
        landmarks

        % The system model used to simulate the platform movement and
        % sensor measurements.
        systemModel;

        % The controller
        platformController;

        % The store of the ground truth state and observations
        xTrueStore;
        timeStore;

    end

    methods(Access = public)
        
        function obj = Simulator(config)
            % Simulator Constructor for Simulator
            %
            % Syntax:
            %   simulator = Simulator(config)
            %
            % Description:
            %   Creates an instance of a Simulator object. The system model
            %   and waypoint controllers are constructed at this time.
            %
            % Inputs:
            %   config - (struct)
            %       The configuration structure
            %
            % Outputs:
            %   simulator - (handle)
            %       An instance of a Simulator

            % Call the base class
            obj@ebe.core.EventBasedSimulator(config);

            % Instantiate the system model
            obj.systemModel = cw1.drivebot.SystemModel(config, config.perturbWithNoise);

            % Instantiate the controller.
            obj.platformController = cw1.drivebot.WaypointController(obj.config);
        end
        
        % Get the ground truth
        function x = xTrue(obj)
            % XTRUE Return the true state of the platform
            %
            % Syntax:
            %   x = simulator.xTrue()
            %
            % Description:
            %   Returns the true state of the platform.
            %
            % Outputs:
            %   x - (3x1 double)
            %       The platform state (x,y,theta)
            x = obj.x;    
        end

        function [timeHistory, xTrueHistory] = history(obj)
            % XTRUE Return the history of the true platform state over time
            %
            % Syntax:
            %   [timeHistory, xTrueHistory] = simulator.history()
            %
            % Description:
            %   Internally, the simulator stores the history of the
            %   platform state over time which can be accessed. This should
            %   be considered deprecated in favour of the results
            %   accumulator.
            %
            % Outputs:
            %   timeHistory - (1xT double)
            %       The simulation time at which each platform state was
            %       recorded.
            %   xTrueHistory - (3xT double)
            %       The platform state (x,y,theta).

            timeHistory = obj.timeStore;
            xTrueHistory = obj.xTrueStore;
        end

        function start(obj)
            % START Start the simulator
            %
            % Syntax:
            %   simulator.start()
            %
            % Description:
            %   Start the simulator. This includes clearing any history
            %   values, starting the platform controller and clearing the
            %   event queue. If the landmark configuration is random, the
            %   landmark positions are also randomly drawn at this time.
            %
            %   The method also kick starts the simulation by scheduling a
            %   call to the initialize handler.

            % Call the base class start functions
            start@ebe.core.EventBasedSimulator(obj);

            obj.platformController.start();

            % Prestore the results
            obj.timeStore = NaN(1, 0);%obj.config.maximumStepNumber + 1);
            obj.xTrueStore = NaN(length(obj.config.x0));%, obj.config.maximumStepNumber + 1);

            % Set up the process noise
            obj.sigmaUSqrtm = diag(obj.config.platform.controller.odomSigma);
            obj.sigmaUSqrtm(3, 3) = deg2rad(obj.sigmaUSqrtm(3,3));
            obj.sigmaU = obj.sigmaUSqrtm * obj.sigmaUSqrtm';

            % Set the initialization callbabk
            obj.eventGeneratorQueue.insert(0, @obj.initialize);

            % If the scenario was not loaded, reject
            if ((isfield(obj.config, 'scenario') == false) || (isstruct(obj.config.scenario) == false))
                return;
            end

            obj.scenario = obj.config.scenario;

            % Are there any landmarks?
            if (isfield(obj.scenario, 'landmarks') == false)
                return
            end

            % Get the landmarks
            slamLandmarks = obj.scenario.landmarks.slam;

            % Just handle random case for now
            if (strcmp(slamLandmarks.configuration, 'random') == true)
                lms = [slamLandmarks.x_min;slamLandmarks.y_min] + ...
                    [slamLandmarks.x_max-slamLandmarks.x_min;slamLandmarks.y_max-slamLandmarks.y_min] .* ...
                    rand(2, slamLandmarks.numLandmarks);
            else
                lms = slamLandmarks.landmarks';
            end

            obj.landmarks = lms;
        end

        function scenario = getScenario(obj)
            % GETSCENARIO Return the scenario used inside the simulator.
            %
            % Syntax:
            %   scenario = simulator.getScenario()
            %
            % Description:
            %   Return the part of the configuration which describes the
            %   scenario (landmarks, GPS occluders, and sensor
            %   configurations).
            %
            % Outputs:
            %   scenario - (struct)
            %       The part of the configuration structure related to
            %       scenarios. (This is a bit of a mess..)

            scenario = obj.scenario;
        end

    end

    methods(Access = protected)

        function handlePredictForwards(obj, dT, ~)
            % HANDLEPREDICTFORWARDS Update the platform position to the
            % next simulation time
            %
            % Syntax:
            %   simulator.handlePredictForwards(dt, DT);
            %
            % Description:
            %   Using the current platform state and controller values,
            %   predict a time dT into the future.
            %
            % Inputs:
            %   dT - (double)
            %       Length of the prediction step.
            %   T - (double)
            %       Simulator wallclock time (not used).

            obj.x = obj.systemModel.predictState(obj.x, obj.u, dT);
        end

        function initialize(obj)
            % INITIALIZE The simulator initialize handler.
            %
            % Syntax:
            %   simulator.initialize();
            %
            % Description:
            %   This is the first step made by the simulator. It
            %   initializes the platform state to its initial conditions
            %   and emits an initialization event. Odomotry is first set.
            %   Initial scheduling times to update odometry and simulate
            %   all enabled sensors are also scheduled.
            %
            %   This handler should be the first thing called by the
            %   simulator.

            assert(obj.stepNumber == 0)

            % Initialize the ground truth state
            P0 = obj.config.P0;
            P0(3,3) = P0(3,3) * (pi / 180)^2;
            obj.x = obj.config.x0 + obj.noiseScale * ebe.utils.psd_sqrtm(P0) * randn(size(obj.config.x0));
            obj.initialized = true;

            % Construct and post the init event
            event = ebe.core.Event(obj.currentTime, 'init', obj.config.x0, obj.config.P0);
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);

            % Force the first odometry event; this also schedules
            % subsequent events
            obj.updateOdometry();

            % Schedule the timeout event
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.heartbeatPeriod, ...
                    @obj.generateHeartbeat);

            % If no scenario is defined, we have no sensors to simulate
            if ((isfield(obj.config, 'scenario') == false) || (isfield(obj.config.scenario, 'sensors') == false))
                return
            end

            % Schedule the bearing measurement
            if (isfield(obj.config.scenario.sensors, 'bearing')  == true)  && (obj.config.scenario.sensors.bearing.enabled == true)
                obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.bearing.measurementPeriod, ...
                    @obj.predictBearingObservations);
            end

            % Schedule the GPS measurement
            if (isfield(obj.config.scenario.sensors, 'gps') == true) && (obj.config.scenario.sensors.gps.enabled == true)
                obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.gps.measurementPeriod, ...
                    @obj.predictGPSObservation);
            end

            % Schedule the compass measurement
            if (isfield(obj.config.scenario.sensors, 'compass') == true) && (obj.config.scenario.sensors.compass.enabled == true)
                obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.compass.measurementPeriod, ...
                    @obj.predictCompassObservation);
            end


             % Schedule the SLAM measurement
            if (isfield(obj.config.scenario.sensors, 'slam')  == true)  && (obj.config.scenario.sensors.slam.enabled == true)
                obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.slam.measurementPeriod, ...
                    @obj.predictSLAMObservations);
            end      
        end

        function updateOdometry(obj)
            % UPDATEODOMETRY Update the platform odometry
            %
            % Syntax:
            %   simulator.updateOdometry();
            %
            % Description:
            %   Call the waypoint controller and work out the new control
            %   inputs for the simulator. Create a noise corrected odom
            %   event which stores this value.
            %
            %   The odometry event data is the noise-corrupted observation
            %   of u (speed, angular velocity) together with the covariance
            %   matrix of the noises.
            %
            %   If the controller has reached the end of the waypoint list
            %   and repeat is not enabled, this method will flag the
            %   simulator that the simulation session is finished.
            %
            %   Outgoing events are emitted and the next call to this
            %   method is scheduled.

            % Get the control input from the controller
            obj.u = obj.platformController.computeControlInputs(obj.x);

            % If this is empty, inform the simulator it's time to stop.
            if (isempty(obj.u))
                obj.carryOnRunning = false;
                return
            end

            % Add noise
            u = obj.u + obj.noiseScale * obj.sigmaUSqrtm * randn(3, 1);

            % Construct and post the event
            event = ebe.core.Event(obj.currentTime, 'odom', u, obj.sigmaU);
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);

            % Schedule next call to this method
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.platform.controller.odomUpdatePeriod, ...
                @obj.updateOdometry);            
        end

        function predictBearingObservations(obj)
            % PREDICTBEARINGOBSERVATIONS Predict the observations from the
            % bearing sensors.
            %
            % Syntax:
            %   simulator.predictBearingObservations();
            %
            % Description:
            %   For each bearing sensor, check if the platform is visible.
            %   If it is, create a bearing observation event. If it is not
            %   detected, a null_obs event is generated.
            %
            %   The bearing observation consists of a list of all bearing
            %   measurements, a covariance matrix (the noise is assumed to
            %   be the same for all sensors) and the list of IDs of bearing
            %   sensors which saw the platform.
            %
            %   Outgoing events are emitted and the next call to this
            %   method is scheduled.

            % Get the list of sensors which can see the platform
            sensorIDs = obj.isDetectedByBearingSensors();

            if (isempty(sensorIDs) == true)
                % If no sensor can see the platform, create a null
                % observation event
                event = ebe.core.Event(obj.currentTime, 'null_obs');
            else
                % Go through each sensor which can detect the platform and
                % compute the noise-corrupted measurement
                nz = numel(sensorIDs);
                z = zeros(1, nz);
                for s = 1 : nz
                    sensor = obj.map.sensors.bearing.sensors(sensorIDs(s));
                    z(s) = obj.systemModel.predictBearingObservation(obj.x, sensor.position, sensor.orientation);
                end
                % Create the bearing observation event.
                event = ebe.core.Event(obj.currentTime, 'bearing', z, ...
                    obj.config.scenario.sensors.bearing.sigmaR, sensorIDs); 
            end

            % Post the event
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);

            % Schedule next call to this method
            obj.eventGeneratorQueue.insert(obj.currentTime + ...
                obj.config.scenario.sensors.bearing.measurementPeriod, @obj.predictBearingObservations);            
        end

        function predictGPSObservation(obj)
            % PREDICTGPSOBSERVATION Predict the GPS observation
            %
            % Syntax:
            %   simulator.predictGPSObservation();
            %
            % Description:
            %   Check if the plaform state lies inside a GPS occluder. If
            %   it does, generate a null_obs event. If it doesn't corrupt a
            %   noise-corrupted GPS observation event.
            %
            %   The GPS observation event consists of data (x,y) and
            %   covariance.
            %
            %   Outgoing events are emitted and the next call to this
            %   method is scheduled.

            % We get the GPS observation only if the robot is not under an
            % occluder
            if (obj.isWithinGPSOccluder() == true)
                event = ebe.core.Event(obj.currentTime, 'null_obs');
            else
                % Generate the observation
		    [z, ~, ~, R] = obj.systemModel.predictGPSObservation(obj.x);

                % Construct the event
                event = ebe.core.Event(obj.currentTime, 'gps', z, R);
            end

            % Post the event
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);

            % Schedule next call to this method
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.gps.measurementPeriod, @obj.predictGPSObservation);
        end

        function predictCompassObservation(obj)
            % PREDICTCOMPASSOBSERVATION Predict the compass observation
            %
            % Syntax:
            %   simulator.predictCompassObservation();
            %
            % Description:
            %   The compass directly measures the heading of the platform.
            %   It is assumed to be always available
            %
            %   The compass observation event consists of data (theta) and
            %   covariance. 
            %
            %   Outgoing events are emitted and the next call to this
            %   method is scheduled.

            % Estimate the compass
            [z, ~, ~, R] = obj.systemModel.predictCompassObservation(obj.x);

            % Construct and post the event
            event = ebe.core.Event(obj.currentTime, 'compass', z, R);
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);

            % Schedule next call to this method
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.compass.measurementPeriod, ...
                @obj.predictCompassObservation);
        end

        function predictSLAMObservations(obj)
            % SLAMOBSERVATIONS Predict the SLAM observations
            % bearing sensors.
            %
            % Syntax:
            %   simulator.predictSLAMObservations();
            %
            % Description:
            %   Check if each landmark in the SLAM map is visible (within a
            %   detection distance) of the platform. For each visible
            %   landmarks, construct a noise-corrupted observation. If no
            %   landmarks are visible, construct a null_obs event.
            %
            %   If Mk landmarks are visible, the data is a (2xMk) array of
            %   (range, bearing) values. The covariance matrix is 2x2 and
            %   is assumed to be the same for all landmarks. The info field
            %   is a (1xMk) vector of the landmark IDs.
            %
            %   Outgoing events are emitted and the next call to this
            %   method is scheduled.

            % Work out the relative distance to all the robots
            dX = obj.landmarks - obj.x(1:2);
            
            % Squared range to each landmark
            R2 = sum(dX.^2,1);
            R = sqrt(R2);

            % Find all the landmarks in range
            landmarkIDs = find(R <= obj.config.scenario.sensors.slam.detectionRange);
            
            if (isempty(landmarkIDs) == true)
                % If no landmarks are visible, create the empty event
                event = ebe.core.Event(obj.currentTime, 'null_obs');
            else
                % For each landmark wourk out the relative range and
                % bearing
                numLandmarks = length(landmarkIDs);
                z = zeros(2, numLandmarks);
                [z(:, 1), ~, ~, ~, R] = obj.systemModel.predictSLAMObservation(obj.x, obj.landmarks(:, landmarkIDs(1)));
                for l  = 2 : numLandmarks
                    z(:, l) = obj.systemModel.predictSLAMObservation(obj.x, obj.landmarks(:, landmarkIDs(l)));
                end
                % Construct the event
                event = ebe.core.Event(obj.currentTime, 'slam', z, R, landmarkIDs);
            end

            % Post the event
            event.eventGeneratorStepNumber = obj.stepNumber;
            obj.outgoingEvents.insert(event);
            
            % Schedule next call to this method
            obj.eventGeneratorQueue.insert(obj.currentTime + obj.config.scenario.sensors.slam.measurementPeriod, @obj.predictSLAMObservations);
        end

        function generateHeartbeat(obj)
            % GENERATEHEARTBEAT Predict the SLAM observations
            % bearing sensors.
            %
            % Syntax:
            %   simulator.generateHeartbeat();
            %
            % Description:
            %   Generates a heartbeat event. This is periodically scheduled
            %   to force the estimator predictor to run. It contains no
            %   data.

            event = ebe.core.Event(obj.currentTime, 'null_obs');
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

    methods(Access = protected)

        function isInside = isWithinGPSOccluder(obj)
            % ISWITHINGPSOCCLUDER Determines if the platform lies within any
            % GPS occluder
            %
            % Syntax:
            %   isInside = simulator.isWithinGPSOccluder();
            %
            % Description:
            %   Check every occluder and return if the platform's (x,y)
            %   position lies within at least one of them.
            %
            % Outputs:
            %   isInside - (bool)
            %       True if the platform lies within a GPS occluder

            % If no occluders are defined, we are never in an occluder
            if (isfield(obj.scenario.sensors.gps, 'occluders') == false)
                isInside = false;
                return
            end

            % Go through all occluders and terminate early if we find an
            % occluder within which the platform lies.
            occluders = obj.scenario.sensors.gps.occluders;

            isInside = false;
            for i = 1:length(occluders)
                occluder = occluders(i);
                if obj.x(1) >= occluder.x_min && obj.x(1) <= occluder.x_max && ...
                   obj.x(2) >= occluder.y_min && obj.x(2) <= occluder.y_max
                    isInside = true;
                    return;
                end
            end
        end
        
        function [sensorIDs] = isDetectedByBearingSensors(obj)
            % ISDETECTEDBYBEARINGSENSORS Determines which, if any, of the
            % bearing sensors can see the platform.
            %
            % Syntax:
            %   isInside = simulator.isWithinGPSOccluder();
            %
            % Description:
            %   Check every bearing sensor and if the platform lies in the
            %   detection cone (distance and angle) place it on the list of
            %   sensorIDs.
            %
            % Outputs:
            %   sensorIDs - (1xN int array)
            %       List of sensor IDs where the platform is visible.

            % Clear the return array
            sensorIDs = [];

            % Get the set of bearing sensors
            sensors = obj.map.sensors.bearing.sensors;

            % Go through each sensor
            for s = 1:length(sensors)
                sensor = sensors(s);
                sensorPos = sensor.position;
                dx = obj.x(1) - sensorPos(1);
                dy = obj.x(2) - sensorPos(2);
                distance = hypot(dx, dy);
                
                % Filter off of detection range and angle.
                if distance <= sensor.detectionRange
                    angleToPoint = atan2(dy, dx);
                    relativeAngle = mod(angleToPoint - deg2rad(sensor.orientation) + pi, 2*pi) - pi;
                    
                    if abs(relativeAngle) <= deg2rad(sensor.detectionAngle)/2
                        sensorIDs(end + 1) = s;
                    end
                end
            end
        end
    end
end
