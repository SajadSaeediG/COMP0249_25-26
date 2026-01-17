classdef SystemModel < handle

    properties(GetAccess = public)

        % Continuous time process model
        Fc;
        Fd;

        % Continuous time process noise model
        Qc;
        Qd;
        QdSqrtm;

        % Observation model
        HGPS;

        % GPS
        RGPS;
        RGPSSqrtm;

        % Bearing
        RBearing;
        RBearingSqrt;

        cachedDT;

    end

    properties(Access = protected)
        config;
        perturbWithNoise;
    end

    methods(Access = public)

        function obj = SystemModel(config, perturbWithNoise)
            obj.config = config;

            if (nargin == 1)
                obj.perturbWithNoise = false;
            else
                obj.perturbWithNoise = perturbWithNoise;
            end
            
            obj.setupModels();

            obj.cachedDT = NaN;
        end

        function [x, F, Q] = predictState(obj, x, dT)

            if (dT ~= obj.cachedDT)
                [obj.Fd, obj.Qd] = ebe.utils.continuousToDiscrete(obj.Fc, ...
                    obj.Qc, dT);
                if (obj.perturbWithNoise == true)
                    obj.QdSqrtm = sqrtm(obj.Qd);
                end
                obj.cachedDT = dT;
            end

            x = obj.Fd * x;

            if (obj.perturbWithNoise == true)
                x = x + obj.QdSqrtm * randn(size(x));
            end

            if (nargout == 3)
                F = obj.Fd;
                Q = obj.Qd;
            end
        end

        function [z, H, R] = predictGPSObservation(obj, x)

           z = obj.HGPS * x;

            if (obj.perturbWithNoise == true)
                z = z + obj.RGPSSqrtm * randn(size(z));
            end

            if (nargout == 3)
                H = obj.HGPS;
                R = obj.RGPS;
            end
        end

        function [z, H, R] = predictBearingObservation(obj, x, sensorXY, sensorTheta)

            dx = x(1) - sensorXY(1);
            dy = x(3) - sensorXY(2);
            z = atan2(dy, dx) - deg2rad(sensorTheta);
            z=atan2(sin(z), cos(z));

            if (obj.perturbWithNoise == true)
                z = z + obj.RBearingSqrt * randn(size(z));
            end

            z = g2o.stuff.normalize_theta(z);

            if (nargout == 3)
                d2 = dx^2 + dy^2;
                H = [-dy/d2 0 dx/d2 0];
                R = obj.RBearing;
            end
        end

    end

    methods(Access = protected)

        function setupModels(obj)

            % Set up the continuous time process model
            obj.Fc=[0 1 0 0;
                -obj.config.alpha -obj.config.beta 0 0;
                0 0 0 1;
                0 0 -obj.config.alpha -obj.config.beta];

            obj.Qc = [0 0 0 0;
                       0 1 0 0;
                       0 0 0 0;
                       0 0 0 1] * obj.config.sigmaQ^2;

            % Set up the observation
            obj.HGPS = [1 0 0 0;0 0 1 0];

            obj.RGPS = eye(2) * obj.config.scenario.sensors.gps.sigmaR^2;

            if (obj.perturbWithNoise == true)
                obj.RGPSSqrtm = obj.config.scenario.sensors.gps.sigmaR;
            end

            if (isfield(obj.config.scenario.sensors, 'bearing'))
                obj.RBearing = deg2rad(obj.config.scenario.sensors.bearing.sigmaR)^2;
                if (obj.perturbWithNoise == true)
                    obj.RBearingSqrt = deg2rad(obj.config.scenario.sensors.bearing.sigmaR);
                end
            end
        end
    end

end