classdef EKFSLAMSystem < cw1.drivebot.SLAMSystemBase
    % EKFSLAMSystem summary of EKFSLAMSystem
    %
    % This class implements an EKF-SLAM system for the drivebot. It is
    % identical to the code in Lab 03, except for the fact that it has now
    % been rewwritten to use the SLAMSystemBase class for common
    % operations.
    %
    % See also:
    %    SLAMSYSTEMBASE

    properties(Access = protected)

        % Although in the lectures we presented three types of states
        % (estimated, predicted, partial) from an implementation point of
        % view it's easier just to store the one set
        x;
        P;

        % The system model
        systemModel;
    end

    methods(Access = public)

        function obj = EKFSLAMSystem(config)
            % SLAMSystem Constructor for SLAMSystem
            %
            % Syntax:
            %   slamSystem = SLAMSystem(config)
            %
            % Description:
            %   Creates an instance of a SLAMSystem object. The system model
            %   and is constructed at this time and a first set of event
            %   handlers are scheduled.
            %
            % Inputs:
            %   config - (struct)
            %       The configuration structure
            %
            % Outputs:
            %   slamSystem - (handle)
            %       An instance of a SLAMSystem

            % Call base class
            obj@cw1.drivebot.SLAMSystemBase(config);

            % Set up the discrete time system for prediction
            obj.systemModel = cw1.drivebot.SystemModel(config);

            % Set the name
            obj.setName('ekf-slam');
        end

        function start(obj)
            % START Start the SLAM system
            %
            % Syntax:
            %   slamSystem.start()
            %
            % Description:
            %   Start the simulator. This includes clearing any results
            %   history and any map.

            start@cw1.drivebot.SLAMSystemBase(obj);

            % Set the dictionary which maps landmark ID to coefficient in
            % the state estimate.
            obj.landmarkIDStateVectorMap = configureDictionary("uint32", "double");

            % Set up initial store of results
            obj.timeStore = [];
            obj.xStore = zeros(cw1.drivebot.SystemModel.NP, 0);
            obj.PStore = zeros(cw1.drivebot.SystemModel.NP, 0);
        end

        function [x,P] = platformEstimate(obj)
            x = obj.x(1:cw1.drivebot.SystemModel.NP);
            P = obj.P(1:cw1.drivebot.SystemModel.NP, 1:cw1.drivebot.SystemModel.NP);
        end
        
        function [T, X, PX] = platformEstimateHistory(obj)
            T = obj.timeStore;
            X = obj.xStore;
            PX = obj.PStore;
        end
        
        function [m, Pmm, landmarkIds] = landmarkEstimates(obj)
            % LANDMARKESTIMATES Return the current mean and covariance of
            % each landmark estimate.
            %
            % Syntax:
            %   [m, Pmm, landmarkIds] = slamSystem.platformEstimate()
            %
            % Description:
            %   Return the current estimates of the landmarks and the IDs.
            %   The landmark covariances are just the blocks on the
            %   diagonals. The full landmark covariance block is not
            %   returned.
            %
            %   At a given time, there are Nk landmarks. The dimension of
            %   each landmark is 2
            %
            % Outputs:
            %   m - (2xN_k vector)
            %       A column vector which contains the estimated mean of
            %       each landmark position.
            %   PX - (2x2xN_k double psd matrix)
            %       A three dimensional matrix which stores the landmark
            %       estimates. The covariance of landmark 4, for example,
            %       is given by P(:,:,4)

            % Get the number of landmarks
            landmarkIds = keys(obj.landmarkIDStateVectorMap);
            numberOfLandmarks = numel(landmarkIds);

            % Build the mean and covariance matrices and extract
            m = NaN(cw1.drivebot.SystemModel.NL, numberOfLandmarks);
            Pmm = NaN(cw1.drivebot.SystemModel.NL, cw1.drivebot.SystemModel.NL, numberOfLandmarks);
            
            for l = 1 : numberOfLandmarks
                landmarkId = landmarkIds(l);
                offset = lookup(obj.landmarkIDStateVectorMap, landmarkId);
                idx = offset + [1;2];
                m(:, l) = obj.x(idx);
                Pmm(:, :, l) = obj.P(idx, idx);
            end
            
        end
    end

    methods(Access = protected)

        function success = handleNoPrediction(obj)
            % HANDLENOPREDICTION Handle the case when no prediction is
            % needed.
            %
            % Syntax:
            %   slamSystem.handleNoPrediction();
            %
            % Description:
            %   This method is called when the time difference between two
            %   events is so small that running the time prediction step
            %   isn't required.

            success = true;
        end

        function success = handleNoUpdate(obj, ~)
            % HANDLENOUPDATE Handle the case when no update is needed.
            %
            % Syntax:
            %   slamSystem.handleNoPrediction();
            %
            % Description:
            %   This method is called when no update to the SLAM system is
            %   required. For example, this is the case with heartbeat and
            %   null_obs events.

            success = true;
        end

        function success = handlePredictForwards(obj, dT)
            % HANDLEPREDICTFORWARDS Predict the system state forwards
            %
            % Syntax:
            %   slamSystem.handlePredictForwards(dT);
            %
            % Description:
            %   Predict the system state forward by a time step dT. The
            %   process model is called with the current platform state
            %   estimate and odometry value. The estimated state and
            %   Jacobians are returned. The covariance estimate is updated
            %   from the Jacobians.
            %
            % Inputs:
            %   dT - (double)
            %       Length of the prediction step.

            tStart = tic;

            % Convenience value
            NP = cw1.drivebot.SystemModel.NP;

            % Update the platform state
            [obj.x(1:NP), gradFx, gradFv] = obj.systemModel.predictState(obj.x(1:NP), obj.u, dT);

            % Update the top left block of the platform state from
            % gradF*Pxx*gradF'+gradV*Q*gradV'
            obj.P(1:NP,1:NP) = gradFx * obj.P(1:NP, 1:NP) * gradFx' + gradFv * obj.covU * gradFv';

            % Do the platform landmark-prediction blocks
            obj.P(1:NP, NP+1:end) = gradFx * obj.P(1:NP, NP+1:end);
            obj.P(NP+1:end, 1:NP) = obj.P(1:NP, NP+1:end)';

            tStop = toc(tStart);

            % Record the time the optimization was called
            obj.performanceData.add('ekf.pr.time', obj.currentTime);
            obj.performanceData.add('ekf.pr.dt', tStop);
            obj.performanceData.add('ekf.pr.str', ...
                numel(obj.x)+numel(obj.x)*(numel(obj.x)-1)*0.5);

            success = true;
        end

        function success = handleInitializationEvent(obj, event)
            % HANDLEINITIALIZATIONEVENT Handle the initialization event.
            %
            % Syntax:
            %   slamSystem.handleInitializationEvent(event);
            %
            % Description:
            %   Handle the initialization event. data is assumed to be x0,
            %   and covariance P0. The initialized flag is also set to true
            %   to enable the estimator to run in full mode.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The intialization event.
            %
            % See Also:
            %   SIMULATOR

            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            obj.x = event.data;
            obj.P = event.covariance;
            obj.initialized = true;
            success = true;
        end

        % Handle a set of measurements of landmarks
        function success = handleSLAMObservationEvent(obj, event)
            % HANDLESLAMOBSERVATIONEVENT Handle the SLAM observation event.
            %
            % Syntax:
            %   slamSystem.handleSLAMObservationEvent(event);
            %
            % Description:
            %   Handle SLAM observations. This follows the code in the
            %   lectures: any known landmarks are updated first. Unknown
            %   landmarks are augmented at the end. Each landmark is
            %   processed individually rather than batching them together.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The SLAM observation event.
            %
            % See Also:
            %   SIMULATOR
            
            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            % Store useful values
            NL = cw1.drivebot.SystemModel.NL;
            NP = cw1.drivebot.SystemModel.NP;

            % Get the list of landmarks we know about
            knownLandmarkIDs = obj.landmarkIDStateVectorMap.keys();

            % Find the intersection of the observed landmarks and the ones
            % known in the map; this produces a list of known landmarks
            % which have been observed this time.
            [existingLandmarks, idx] = intersect(event.info, knownLandmarkIDs);

            obj.performanceData.add('ekf.slm.time', obj.currentTime);
            obj.performanceData.add('ekf.slm.nel', numel(idx));

            tStart = tic;

            % Update all the known landmarks
            for o = 1 : numel(idx)
                % Look up the landmark and figure out its position
                offset = lookup(obj.landmarkIDStateVectorMap, existingLandmarks(o));
                landmarkIdx = offset + (1:NL);

                % Predicted observation and Jacobians
                [zPred, gradHx, gradHm, gradHw] = ...
                    obj.systemModel.predictSLAMObservation(obj.x(1:NP), ...
                    obj.x(landmarkIdx));

                % Work out the innovation, including angle wrapping
                nu = event.data(:, idx(o)) - zPred;
                nu(2) = atan2(sin(nu(2)), cos(nu(2)));

                % Assemble the observation matrix
                HS = zeros(2, numel(obj.x));
                HS(:, 1:NP) = gradHx;
                HS(:, landmarkIdx) = gradHm;

                % Kalman filter update steps
                C = obj.P * HS';
                S = HS * C + gradHw * event.covariance() * gradHw';
                K = C / S;
                obj.x = obj.x + K * nu;
                obj.P = obj.P - K * S * K';

                % Wrap the heading estimate
                obj.x(3) = atan2(sin(obj.x(3)), cos(obj.x(3)));

            end

            obj.performanceData.add('ekf.slm.el_dt', toc(tStart));

            % Find the mutual complement  of the observed landmarks and the
            % ones known in the map; this produces a list of new landmarks
            % which have not been seen before.
            [newLandmarks, idx] = setdiff(event.info, existingLandmarks);

            obj.performanceData.add('ekf.slm.nnl', numel(newLandmarks));

            tStart = tic;

            % Augment all the known landmarks
            for o = 1 : numel(newLandmarks)

                % Figure out the index for the new state estimate and
                % insert it into the look up table
                stateDim = length(obj.x);
                landmarkIdx = stateDim + (1:NL);
                obj.landmarkIDStateVectorMap = ...
                    insert(obj.landmarkIDStateVectorMap, newLandmarks(o), stateDim);

                % Use the inverse observation model to estimate the
                % landmark position and compute the Jacobian
                [m, gradGx, gradGw] = ...
                    obj.systemModel.predictLandmarkFromSLAMObservation(obj.x, event.data(:, idx(o)));

                % Expand the state with the initial estimate of the
                % landmark position
                obj.x(landmarkIdx) = m;

                % Add the predicted landmark covariance. We do this first
                % because it forces obj.P to be resized to the correct
                % dimension
                obj.P(landmarkIdx, landmarkIdx) = gradGx * obj.P(1:NP, 1:NP) * gradGx' ...
                    + gradGw * event.covariance() * gradGw';

                % Compute the cross correlation
                obj.P(landmarkIdx, 1:end-NL) = gradGx * obj.P(1:NP, 1:end-NL);

                % Copy over tranposed version
                obj.P(1:end-NL,landmarkIdx) = obj.P(landmarkIdx, 1:end-NL)';
            end

            obj.performanceData.add('ekf.slm.nl_dt', toc(tStart));

            success = true;
        end

        function success = handleGPSObservationEvent(obj, event)
            % HANDLEGPSOBSERVATIONEVENT Handle the GPS observation event.
            %
            % Syntax:
            %   slamSystem.handleGPSObservationEvent(event);
            %
            % Description:
            %   Handle the GPS observation.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The GPS observation event.
            %
            % See Also:
            %   SIMULATOR

            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            % Predicted observation and Jacobians
            [zPred, gradHx, gradHw] = ...
                obj.systemModel.predictGPSObservation(obj.x(1:2));

            % Compute the innovation
            nu = event.data - zPred;

            % Assemble the observation matrix
            HS = zeros(2, numel(obj.x));
            HS(:, 1:cw1.drivebot.SystemModel.NP) = gradHx;

            % Kalman Filter Update
            C = obj.P * HS';
            S = HS * C + gradHw * event.covariance() * gradHw';
            W = C / S;
            obj.x = obj.x + W * nu;
            obj.P = obj.P - W * S * W';

            % Wrap the heading estimate
            obj.x(3) = atan2(sin(obj.x(3)), cos(obj.x(3)));

            success = true;
        end

        function success = handleCompassObservationEvent(obj, event)
            % HANDLECOMPASSOBSERVATIONEVENT Handle the compass observation event.
            %
            % Syntax:
            %   slamSystem.handleCompassObservationEvent(event);
            %
            % Description:
            %   Handle the GPS observation.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The compass observation event.
            %
            % See Also:
            %   SIMULATOR

            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            % Predicted observation and Jacobians
            [zPred, gradHx, gradHw] = ...
                obj.systemModel.predictCompassObservation(obj.x(1:cw1.drivebot.SystemModel.NP));

            % Compute the innovation
            nu = event.data - zPred;
            nu = atan2(sin(nu), cos(nu));

            % Expand to the full state
            HS = zeros(1, numel(obj.x));
            HS(:, 1:cw1.drivebot.SystemModel.NP) = gradHx;

            % Kalman Filter Update
            C = obj.P * HS';
            S = HS * C + gradHw * event.covariance() * gradHw';
            W = C / S;
            obj.x = obj.x + W * nu;
            obj.P = obj.P - W * S * W';

            % Wrap the heading estimate
            obj.x(3) = atan2(sin(obj.x(3)), cos(obj.x(3)));

            success = true;
        end


        function success = handleBearingObservationEvent(obj, event)
            % HANDLEBEARINGOBSERVATIONEVENT Handle the bearing observation event.
            %
            % Syntax:
            %   slamSystem.handleBearingObservationEvent(event);
            %
            % Description:
            %   Handle the bearing observation event. For each observation,
            %   we extract the sensor ID and from that work out the bearing
            %   sensor position and orientation. This is then used to carry
            %   out the update.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The bearing observation event.
            %
            % See Also:
            %   SIMULATOR

            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            % Update each measurement separately
            for s = 1 : numel(event.info)

                % Predicted observation and Jacobians
                sensor = obj.map.sensors.bearing.sensors(event.info(s));
                [zPred, gradHx, gradHw] = ...
                obj.systemModel.predictBearingObservation(obj.x(1:2), ...
                    sensor.position, sensor.orientation);

                % Compute the innovation
                nu = event.data(s) - zPred;
                nu = atan2(sin(nu), cos(nu));

                % Expand to full state
                HS = zeros(1, numel(obj.x));
                HS(:, 1:l2s.drivebot.SystemModel.NP) = gradHx;
                
                % Kalman Filter Update
                C = obj.P * HS';
                S = HS * C + gradHw * event.covariance() * gradHw';
                W = C / S;
                obj.x = obj.x + W * nu;
                obj.P = obj.P - W * S * W';

                % Wrap the heading estimate
                obj.x(3) = atan2(sin(obj.x(3)), cos(obj.x(3)));
            end

            success = true;
        end

        function storeStepResults(obj)
            % Store the estimate for the future
            obj.timeStore(:, obj.stepNumber + 1) = obj.currentTime;
            obj.xStore(:, obj.stepNumber + 1) = obj.x(1:cw1.drivebot.SystemModel.NP);
            obj.PStore(:, obj.stepNumber + 1) = diag(obj.P(1:cw1.drivebot.SystemModel.NP, ...
                1:cw1.drivebot.SystemModel.NP));
        end
    end
end