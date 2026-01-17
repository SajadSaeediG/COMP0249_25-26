classdef KalmanFilter < ebe.localization.LocalizationSystem

    properties(Access = protected)

        % Kalman filter prediction
        xPred;
        PPred;

        % Kalman filter estimate
        xEst;
        PEst;

        % The linear system used to predict the state and the observation
        systemModel;

        % The map
        scenario

        % Store of the mean and covariance values
        timeStore;
        xStore;
        PStore;

    end

    methods(Access = public)

        function obj = KalmanFilter(config)

            % Call base class
            obj@ebe.localization.LocalizationSystem(config);

            % Set up the discrete time system for prediction
            obj.systemModel = l1s.pointbot.SystemModel(config);

            % Set up the event handlers
            obj.registerEventHandler('init', @obj.handleInitializationEvent);
            obj.registerEventHandler('null_obs', @obj.handleNoObservation);
            obj.registerEventHandler('heartbeat', @obj.handleNoObservation);
            obj.registerEventHandler('gps', @obj.handleGPSObservationEvent);
            obj.registerEventHandler('bearing', @obj.handleBearingObservationEvent);

            % Set the name
            obj.setName('KalmanFilter');
        end

        function success = start(obj)
            start@ebe.core.EventBasedEstimator(obj);
            obj.timeStore = [];
            obj.xStore = zeros(4, 0);

            % Get the map data
            if (isfield(obj.config, 'scenario'))
                obj.scenario = obj.config.scenario;
            end

            success = true;
        end

        function [T, X, PX] = estimateHistory(obj)
            T = obj.timeStore;
            X = obj.xStore;
            PX = obj.PStore;
        end

        function [x, P] = computeXP(obj)
            x = obj.xEst;
            P = obj.PEst;
        end

    end

    methods(Access = protected)

        function success = handleNoPrediction(obj)
            obj.xPred = obj.xEst;
            obj.PPred = obj.PEst;
            success = true;
        end

        function success = handleNoObservation(obj, event)
            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            fprintf('handleNoObservation(%f): %s\n', obj.currentTime, event.type)

            obj.xEst = obj.xPred;
            obj.PEst = obj.PPred;
            success = true;
        end

        function success = handlePredictForwards(obj, dT)

            %fprintf('handlePredictForwards(%f): dT=%f\n', obj.currentTime, dT)

            [obj.xPred, F, Q] = obj.systemModel.predictState(obj.xEst, dT);

            obj.PPred = F * obj.PEst * F' + Q;

            success = true;
        end

        function success = handleInitializationEvent(obj, event)
            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            %fprintf('handleInitializationEvent(%f)\n', obj.currentTime)

            obj.xEst = event.data;
            obj.PEst = event.covariance;
            obj.initialized = true;
            success = true;
        end

        function success = handleGPSObservationEvent(obj, event)
            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            %fprintf('handleGPSObservationEvent(%f): complete implementation\n', obj.currentTime)

            [zPred, H, R] = obj.systemModel.predictGPSObservation(obj.xPred);

            % Put code down here to complete the implementation of the KF
            nu = event.data - zPred;
            C = obj.PPred * H';
            S = H * C + R;
            W = C / S;
            obj.xEst = obj.xPred + W * nu;
            obj.PEst = obj.PPred - W * S * W';
            success = true;

        end

        function success = handleBearingObservationEvent(obj, event)
            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            %fprintf('handleBearingObservationEvent(%f)\n', obj.currentTime)

            x = obj.xPred;
            P = obj.PPred;

            % Update each measurement separately
            for s = 1 : numel(event.info)
                sensor = obj.scenario.sensors.bearing.sensors(event.info(s));
                [zPred, H, R] = obj.systemModel.predictBearingObservation(x, sensor.position, sensor.orientation);

                % Put code down here to complete the implementation of the KF
                nu = event.data(s) - zPred;
                nu = atan2(sin(nu), cos(nu));
                C = P * H';
                S = H * C + R;
                W = C / S;
                x = x + W * nu;
                P = P - W * S * W';
            end

            obj.xEst = x;
            obj.PEst = P;

            success = true;
        end

        function storeStepResults(obj)
            % Store the estimate for the future
            obj.timeStore(:, obj.stepNumber + 1) = obj.currentTime;
            obj.xStore(:, obj.stepNumber + 1) = obj.xEst;
            obj.PStore(:, obj.stepNumber + 1) = diag(obj.PEst);

        end
    end
end