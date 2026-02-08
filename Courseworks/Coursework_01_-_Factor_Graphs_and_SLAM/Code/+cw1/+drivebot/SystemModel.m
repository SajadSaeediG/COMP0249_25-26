classdef SystemModel < handle
    % SystemModel summary of SystemModel
    %
    % This class stores all of the low-level mathematical functions, such
    % as the process model, observation models and inverse observation
    % models. These equations are presented outside of the context of the
    % specific system (e.g., for a SLAM system, they do not know how many
    % landmarks are in the system, index of those states, etc.)

    properties(Access = public, Constant)
        % Platform state dimension
        NP = 3;
        
        % Landmark dimension
        NL = 2;
    end

    properties(Access = protected)
        % The GPS observation noise covariance matrix
        RGPS;

        % Square root of the GPS observation noise covariance matrix.
        % This is only used by the simulator for sampling the observation
        % noise.
        RGPSSqrtm;

        % Compass covariance
        RCompass;

        % Square root of compass observation noise covariance matrix.
        % This is only used by the simulator for sampling the observation
        % noise.
        RCompassSqrt;

        % Bearing observation noise covariance matrix
        RBearing;

        % Square root of the bearing observation noise covariance matrix.
        % This is only used by the simulator for sampling the observation
        % noise.
        RBearingSqrt;

        % SLAM observation covariance matrix
        RSLAM;

        % Square root of SLAM observation noise covariance matrix.
        % This is only used by the simulator for sampling the observation
        % noise.
        RSLAMSqrt;

        % The system configuration.
        config;

        % A flag to show if noises should be sampled.
        % This flag is only set to true in the simulator.
        perturbWithNoise;
    end

    methods(Access = public)

        function obj = SystemModel(config, perturbWithNoise)
            % SystemModel Constructor for SystemModel
            %
            % Syntax:
            %   obj = SystemModel(config)
            %   obj = SystemModel(config, perturbWithNoise)
            %
            % Description:
            %   Creates an instance of a SystemModel object.
            %
            % Inputs:
            %   config - (struct)
            %       The configuration structure
            %   perturbWithNoise - (optional, bool)
            %       Set to true if noises should be sampled [Default: false]
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a SystemModel
      
            obj.config = config;

            % Handle argument defaults
            if (nargin == 1)
                obj.perturbWithNoise = false;
            else
                obj.perturbWithNoise = perturbWithNoise;
            end
            
            % Do any setup of models
            obj.setupModels();
        end

        function [x, gradFx, gradFv] = predictState(obj, x, u, dT)
            % PREDICTSTATE predict the future state of the platform
            %
            % Syntax:
            %   x = obj.predictState(x, u, dT)
            %   [x, gradFx, gradFv] = obj.predictState(x, u, dT)
            %
            % Description:
            %   Predicts the next state of the platform given the current
            %   state, control input and prediction interval. If three
            %   output arguments are specified, the Jacobians with respect
            %   to platform state and process noise are computed. Note that
            %   even if perturbWithNoise is set to true, it is not applied
            %   in the prediction step. The reason is that the noise is
            %   added to the odometry measurement which is handled
            %   elsewhere.
            %
            % Inputs:
            %   x - (3x1 double)
            %       Platform state (x,y,theta)
            %   u - (double)
            %       The control input (s, omega)
            %   dT - (double)
            %       Length of the prediction step
            %
            % Outputs:
            %   x - (3x1 double)
            %       The predicted platform state (x,y,theta) dT seconds
            %       later.
            %   gradFx - (3x3 double, optional)
            %       The Jacobian of the process model wrt to the
            %       platform state.
            %   gradFv - (3x2 double, optional)
            %       The Jacobian of the process model wrt to the
            %       process noise.

            % Speed times distance
            sDT = u(1) * dT;

            % The M matrix
            M = [cos(x(3)) -sin(x(3)) 0;sin(x(3)) cos(x(3)) 0;0 0 1];

            % Predict forwards using the model
            x = x + dT * M * u;

            % Wrap the angle
            x(3) = atan2(sin(x(3)), cos(x(3)));

            % If requested, compute the Jacobian matrices
            if (nargout == 3)
                gradFx = eye(cw1.drivebot.SystemModel.NP);
                gradFx(1, 3) = -sDT * sin(x(3));
                gradFx(2, 3) = sDT * cos(x(3));
                gradFv = dT * M;
            end
        end

        function [z, gradHx, gradHm, gradHw, R] = predictSLAMObservation(obj, x, mXY)
            % PREDICTSLAMOBSERVATION predict the observation from the SLAM sensor
            %
            % Syntax:
            %   z = obj.predictSLAMObservation(x, mXY)
            %   [z, gradHx, gradHm, gradHw, R] = obj.predictSLAMObservation(x, mXY)
            %
            % Description:
            %   Predict the output of the SLAM sensor (range, bearing)
            %   given that the platform state is x and the landmark state
            %   is mXY. If five output arguments are specified, the
            %   Jacobians of the observation model with respect to the
            %   platform state, landmark state and observation noise,
            %   together with the observation covariance matrix, are
            %   provided.
            %
            % Inputs:
            %   x - (3x1 double)
            %       Platform state (x,y,theta)
            %   mXY - (2x1 double)
            %       Landmark (x,y) position
            %
            % Outputs:
            %   z - (2x1 double)
            %       The predicted SLAM observation (r, beta)
            %   gradHx - (2x3 double, optional)
            %       The Jacobian of the observation model wrt to the
            %       platform state.
            %   gradHm - (2x2 double, optional)
            %       The Jacobian of the observation model wrt to the 
            %       landmark estimate.
            %   gradHw - (2x2 double, optional)
            %       The Jacobian of the observation model wrt to the 
            %       observation noise.
            %   R - (2x2 double, optional)
            %       The covariance of the SLAM observation

            % Work out the relative distance
            dXY = mXY - x(1:2);

            % Range and squared range
            r2 = sum(dXY.^2);
            r = sqrt(r2);

            % Compute noiseless observation
            z = [r;
                atan2(dXY(2), dXY(1)) - x(3)];

            % Add noise if required
            if (obj.perturbWithNoise == true)
                z = z + obj.RSLAMSqrt * randn(size(z));
            end

            % Clamp bearing measurement
            z(2) = atan2(sin(z(2)), cos(z(2)));

            if (nargout >= 3)
                gradHx = [-dXY(1)/r -dXY(2)/r 0;
                    dXY(2)/r2 -dXY(1)/r2 -1];
                gradHm = -gradHx(1:2, 1:2);
                gradHw = eye(2);
                R = obj.RSLAM;
            end
        end

        function [mXY, gradGx, gradGw, R] = predictLandmarkFromSLAMObservation(obj, x, z)
            % PREDICTLANDMARKFROMOBSERVATION the inverse observation model
            % which predicts the landmark state from the platform state and
            % observation.
            %
            % Syntax:
            %   mXY = obj.predictLandmarkFromSLAMObservation(x, z)
            %   [mXY, gradGx, gradGw] = obj.predictLandmarkFromSLAMObservation(x, z)
            %
            % Description:
            %   Predict the output of the SLAM sensor (range, bearing)
            %   given that the platform state is x and the landmark state
            %   is mXY. If four output arguments are specified, the
            %   Jacobians of the inverse observation model with respect to
            %   the platform state and the observation noise, together with the
            %   observation covariance matrix, are computed.
            %
            % Inputs:
            %   x - (3x1 double)
            %       Platform state (x,y,theta)
            %   z - (2x1 double)
            %       The SLAM observation (r, beta)
            %
            % Outputs:
            %   mXY - (2x1 double)
            %       Landmark (x,y) position
            %   gradGx - (2x3 double, optional)
            %       The Jacobian of the observation model wrt to the
            %       platform state.
            %   gradGw - (2x2 double, optional)
            %       The Jacobian of the process model wrt to the
            %       observation noise.
            %   R - (2x2 double, optional)
            %       The covariance of the SLAM observation

            % Angle in world coordinates
            phi = x(3) + z(2);

            % Work out Polar to Carteisan conversion
            mXY = x(1:2) + z(1) * [cos(phi); sin(phi)];

            % Set up the optional outputs if requested
            if (nargout >= 3)
                gradGx = eye(2, 3);
                gradGx(1, 3) = -z(1) * sin(phi);
                gradGx(2, 3) = z(1) * cos(phi);
    
                gradGw = zeros(2, 2);
                gradGw(1, 1) = cos(phi);
                gradGw(1, 2) = -z(1) * sin(phi);
                gradGw(2, 1) = sin(phi);
                gradGw(2, 2) = z(1) * cos(phi);
                R = obj.RSLAM;
            end
        end

        function [z, gradHx, gradHw, R] = predictGPSObservation(obj, x)
            % PREDICTGPSOBSERVATION predict the observation from an idealised GPS sensor
            %
            % Syntax:
            %   z = obj.predictGPSObservation(x)
            %   [x, gradHx, gradHw, R] = obj.predictGPSObservation(x)
            %
            % Description:
            %   Predict the output of the GPS sensor (x, y) from the
            %   platform state (x,y,theta). If four output arguments are
            %   specified, the gradient of the observation model with
            %   respect to the platform state and the observation noise,
            %   together with the observation covariance matrix, are
            %   computed.
            %
            % Inputs:
            %   x - (3x1 double)
            %       Platform state (x,y,theta)
            %
            % Outputs:
            %   z - (2x1 double)
            %       The predicted SLAM observation (r, beta)
            %   gradHx - (2x3 double, optional)
            %       The Jacobian of the observation model wrt to the
            %       platform state.
            %   gradHw - (2x2 double, optional)
            %       The Jacobian of the process model wrt to the 
            %       observation noise.
            %   R - (2x2 double, optional)
            %       The covariance of the SLAM observation

            % Extract the platform position estimate
            z = x(1:2);

            % Add observation noise if requested
            if (obj.perturbWithNoise == true)
                z = z + obj.RGPSSqrtm * randn(size(z));
            end

            % Generate output matrices if requested
            if (nargout >= 2)
                gradHx = eye(2, 3);
                gradHw = eye(2);
                R = obj.RGPS;
            end
        end

        function [z, gradHx, gradHw, R] = predictCompassObservation(obj, x)
            % PREDICTGPSOBSERVATION predict the observation from an idealised compass
            %
            % Syntax:
            %   z = obj.predictCompassObservation(x)
            %   [x, gradHx, gradHw, R] = obj.predictCompassObservation(x)
            %
            % Description:
            %   Predict the output of a compass (theta) from the
            %   platform state (x,y,theta). If four output arguments are
            %   specified, the gradient of the observation model with
            %   respect to the platform state and the observation noise,
            %   together with the observation covariance matrix, are
            %   computed.
            %
            % Inputs:
            %   x - (3x1 double)
            %       Platform state (x,y,theta)
            %
            % Outputs:
            %   z - (double)
            %       The predicted compass heading (theta)
            %   gradHx - (1x3 double, optional)
            %       The Jacobian of the observation model wrt to the
            %       platform state.
            %   gradHw - (2x2 double, optional)
            %       The Jacobian of the process model wrt to the 
            %       observation noise.
            %   R - (2x2 double, optional)
            %       The covariance of the SLAM observation

            % Extract the platform position estimate
            z = x(3);

            % Add observation noise if requested
            if (obj.perturbWithNoise == true)
                z = z + obj.RCompassSqrt * randn;
            end

            % Angle wrap
            z = atan2(sin(z), cos(z));

            % Generate output matrices if requested
            if (nargout >= 2)
                gradHx = [0 0 1];
                gradHw = 1;
                R = obj.RCompass;
            end
        end


        function [z, gradHx, gradHw, R] = predictBearingObservation(obj, x, sensorXY, sensorTheta)
            % PREDICTBEARINGOBSERVATION predict the observation from an bearing sensor
            %
            % Syntax:
            %   z = obj.predictBearingObservation(x, sensorXY, sensorTheta)
            %   [x, gradHx, gradHw, R] = obj.predictBearingObservation(, sensorXY, sensorTheta)
            %
            % Description:
            %   Predict the output of a bearing (delta) from a sensor at a
            %   fixed position and orientation (sensorXY, sensorTheta). If
            %   four output arguments are specified, the gradient of the
            %   observation model with respect to the platform state, the
            %   observation noise are computed together with the
            %   observation covariance matrix.
            %
            % Inputs:
            %   x - (3x1 double)
            %       Platform state (x,y,theta)
            %   sensorXY - (2x1 double)
            %       The sensor (x,y) position
            %   sensorTheta - (double)
            %       The sensor orientation
            %
            % Outputs:
            %   z - (double)
            %       The predicted bearing observation (delta)
            %   gradHx - (1x3 double, optional)
            %       The Jacobian of the observation model wrt to the platform state.
            %   gradHw - (double, optional)
            %       The Jacobian of the observation model wrt to the observation noise.
            %   R - (double, optional)
            %       The covariance of the bearing observation
            dx = x(1) - sensorXY(1);
            dy = x(2) - sensorXY(2);
            z = atan2(dy, dx) - deg2rad(sensorTheta);

            if (obj.perturbWithNoise == true)
                z = z + obj.RBearingSqrt * randn(size(z));
            end

            z = atan2(sin(z), cos(z));

            if (nargout >= 2)
                d2 = dx^2 + dy^2;
                gradHx = [-dy/d2 dx/d2 0];
                gradHw = 1;
                R = obj.RBearing;
            end
        end
    end

    methods(Access = protected)

        function setupModels(obj)
            % Set up the models
            % This method creates elements as needed

            % Validate the configuration file
            assert(isfield(obj.config, 'platform') == true, 'ls3systemmodel:no:platform', ...
                'The configuration is missing the platform configuration');

            if (isfield(obj.config, 'scenario') == false)
                return
            end

            if (isfield(obj.config.scenario, 'sensors') == false)
                return
            end

            % Construct the GPS observation covariance matrix
            if (isfield(obj.config.scenario.sensors, 'gps'))
                obj.RGPS = eye(2) * obj.config.scenario.sensors.gps.sigmaR^2;
    
                if (obj.perturbWithNoise == true)
                    obj.RGPSSqrtm = eye(2) * obj.config.scenario.sensors.gps.sigmaR;
                end
            end

            % Construct the GPS observation covariance matrix
            if (isfield(obj.config.scenario.sensors, 'compass'))
                obj.RCompass = deg2rad(obj.config.scenario.sensors.compass.sigmaR)^2;
    
                if (obj.perturbWithNoise == true)
                    obj.RCompassSqrt = deg2rad(obj.config.scenario.sensors.compass.sigmaR);
                end
            end


            % Construct the SLAM observation covariance matrix
            if (isfield(obj.config.scenario.sensors, 'slam'))
                obj.RSLAM = diag([obj.config.scenario.sensors.slam.sigmaR(1)^2; ...
                    deg2rad(obj.config.scenario.sensors.slam.sigmaR(2))^2]);
                if (obj.perturbWithNoise == true)
                    obj.RSLAMSqrt = ebe.utils.psd_sqrtm_chol(obj.RSLAM);
                end
            end

            % Construct the bearing observation covariance matrix
            if (isfield(obj.config.scenario.sensors, 'bearing'))
                obj.RBearing = deg2rad(obj.config.scenario.sensors.bearing.sigmaR)^2;
                if (obj.perturbWithNoise == true)
                    obj.RBearingSqrt = deg2rad(obj.config.scenario.sensors.bearing.sigmaR);
                end
            end

        end
    end

end