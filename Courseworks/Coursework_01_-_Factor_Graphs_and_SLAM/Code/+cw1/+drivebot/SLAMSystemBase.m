classdef SLAMSystemBase < ebe.slam.SLAMSystem
    % SLAMSystemBase summary of SLAMSystemBase
    %
    % The SLAMSystemBase implements the base for a SLAM system using a
    % Kalman filter or a factor graph. We use a base class to ensure a
    % common structure

    properties(Access = protected)

        % The most recent version of the odometry together with its
        % measurement covariance
        u;
        covU;

        % Map stores landmark ID with the indices in the state vector
        landmarkIDStateVectorMap;

        % The map
        scenario;

        % Time store; we use this to store different times

        % Store of the mean and covariance values
        timeStore;
        xStore;
        PStore;

        % Performance data
        performanceData;
    end

    methods(Access = public)

        function obj = SLAMSystemBase(config)
            % SLAMSystemBase Constructor for SLAMSystemBase
            %
            % Syntax:
            %   slamSystemBase = SLAMSystemBase(config)
            %
            % Description:
            %   Creates an instance of a SLAMSystemBase object. This is the
            %   base class for different SLAM implementations. This class
            %   is abstract and cannot be directly instantiated.
            %
            % Inputs:
            %   config - (struct)
            %       The configuration structure
            %
            % Outputs:
            %   slamSystemBase - (handle)
            %       An instance of a SLAMSystemBase

            % Call base class
            obj@ebe.slam.SLAMSystem(config);

            % Set the name
            obj.setName('SLAMSystemBase');
        end

        function start(obj)
            % START Start the SLAM system
            %
            % Syntax:
            %   slamSystem.start()
            %
            % Description:
            %   Start the simulator. This includes clearing any results
            %   history and any map. The event handlers are also
            %   registered.

            start@ebe.slam.SLAMSystem(obj);

            % Set up the event handlers
            obj.registerEventHandler('init', @obj.handleInitializationEvent);
            obj.registerEventHandler('null_obs', @obj.handleNoUpdate);
            obj.registerEventHandler('gps', @obj.handleGPSObservationEvent);
            obj.registerEventHandler('slam', @obj.handleSLAMObservationEvent);
            obj.registerEventHandler('odom', @obj.handleUpdateOdometryEvent);
            obj.registerEventHandler('compass', @obj.handleCompassObservationEvent);

            % Get the map data
            if (isfield(obj.config, 'scenario'))
                obj.scenario = obj.config.scenario;
            end

            % Set up the performance data
            obj.performanceData = ebe.utils.PerformanceData();
        end

        function performanceData = getPerformanceData(obj)
            % GETPERFORMANCEDATA Return the performance data structure.
            %
            % Syntax:
            %   performanceData = slamSystem.getPerformanceData();
            %
            % Return:
            %   performanceData - (ebe.utils.PerformanceData)
            %       The performance data stored with this object.
            %
            % See also:
            %   EBE.UTILS.PERFORMANCEDATA

            performanceData = obj.performanceData;
        end
    end

    methods(Access = protected, Sealed)

        function success = handleUpdateOdometryEvent(obj, event)
            % HANDLEUPDATEODOMETRYEVENT Handle receving new odometry
            % information.
            %
            % Syntax:
            %   slamSystem.handleUpdateOdometryEvent(event);
            %
            % Description:
            %   Handle the odometry event. The internal value of odometry
            %   is updated to the value in the event. The odometry is held
            %   constant until the next handle event is received.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The odometry event.
            %
            % See Also:
            %   SIMULATOR

                assert(obj.stepNumber == event.eventGeneratorStepNumber);
                obj.u = event.data;
                obj.covU = event.covariance;
                success = true;
        end

        function storePerformanceData(obj, operation, value)
            % STOREPERFORMANCEDATA Store performance data for subsequent analysis.
            %
            % Syntax:
            %   slamSystem.storePerformanceData(operation, value);
            %
            % Description:
            %   For various questions you will need to store data such as
            %   the time required and size of matrices. This method
            %   supports a key-value pairing to do this
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The odometry event.

        end

    end

    methods(Access = protected, Abstract)

        success = handleNoPrediction(obj);
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

         success = handleNoUpdate(obj, event);
            % HANDLENOUPDATE Handle the case when no update is needed.
            %
            % Syntax:
            %   slamSystem.handleNoPrediction();
            %
            % Description:
            %   This method is called when no update to the SLAM system is
            %   required. For example, this is the case with heartbeat and
            %   null_obs events.

        success = handlePredictForwards(obj, dT);
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

            % Convenience value

        success = handleInitializationEvent(obj, event);
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

        % Handle a set of measurements of landmarks
        success = handleSLAMObservationEvent(obj, event);
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

         success = handleGPSObservationEvent(obj, event);
            % HANDLEGPSOBSERVATIONEVENE Handle the GPS observation event.
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

         success = handleCompassObservationEvent(obj, event);
            % HANDLECOMPASSOBSERVATION Handle the compass observation event.
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

         success = handleBearingObservationEvent(obj, event);
            % HANDLECBEARINGOBSERVATIONEVENT Handle the bearing observation event.
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
    end
end