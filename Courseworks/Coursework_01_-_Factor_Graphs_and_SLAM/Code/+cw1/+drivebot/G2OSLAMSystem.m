classdef G2OSLAMSystem < cw1.drivebot.SLAMSystemBase
    % G2OSLAMSystem summary of G2OSLAMSystem
    %
    % This class implements the SLAM system using the G2O library.
    %
    % See also:
    %    SLAMSYSTEMBASE

    properties(Access = public, Constant)
        % Platform state dimension
        NP = 3;
        
        % Landmark dimension
        NL = 2;
        
        % Initial cache size; might help a bit with performance
        INITIAL_CACHE_SIZE = 10000;
    end

    properties(Access = public)

        % The most recently created vehicle vertex.
        currentPlatformVertex;
        
        % The set of all vertices associated with the vehicle state over
        % time.
        platformVertices;
        platformVertexId;
        
        % The set of all prediction edges. These are removed from the graph
        % afterwards if we don't use prediction
        processModelEdges;
        numProcessModelEdges;

        % Flag to run the detailed graph validation checks
        validateGraphOnInitialization;
        
        % The graph used for performing estimation.
        graph;
        
        % The optimization algorithm
        optimizationAlgorithm;

        % Settings to limit the number of landmark observtions; used in Q3a
        limitObservationsPerLandmark;
        maxObservationsPerLandmark;

        % Setting to fix older vehicle poses; used in Q3b.
        fixOlderPlatformVertices;
        unfixedTimeWindow;

    end

    methods(Access = public)

        function obj = G2OSLAMSystem(config)
            % SLAMSystem Constructor for SLAMSystem
            %
            % Syntax:
            %   slamSystem = G2OSLAMSystem(config)
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

            % Disable validation checks for speed
            obj.validateGraphOnInitialization = false;

            % Initially disable graph pruning
            obj.limitObservationsPerLandmark = false;
            obj.maxObservationsPerLandmark = NaN;

            % Disable older vertex fixing
            obj.fixOlderPlatformVertices = false;
            obj.unfixedTimeWindow = NaN;

            % Set the name
            obj.setName('g2o-slam');
        end

        function setMaxObservationsPerLandmark(obj, maxObservationsPerLandmark)
            % SETMAXOBSERVATIONSPERLANDMARK Set the maximum number of
            % observations the graph should maintain per landmark
            %
            % Syntax:
            %   slamSystem.setMaxObservationsPerLandmark(maxObservationsPerLandmark)
            %
            % Description:
            %   The basic graph continues to add edges to the landmark,
            %   which causes computational and storage costs to keep
            %   increasing. This method specifies an upper limit. If it is
            %   exceeded, landmarks are randomly removed. Note the first
            %   and most recent observations are always kept.
            %
            % Inputs:
            %   maxObservationsPerLandmark - (int)
            %       The maximum number of obsevations per landmark

            obj.maxObservationsPerLandmark = maxObservationsPerLandmark;
            obj.limitObservationsPerLandmark = true;
        end

        function setFixOlderPlatformVertices(obj, unfixedTimeWindow)
            % SETFIXOLDERPLATFORMVERTICES Set the maximum time a platform vertex remains unfixed.
            %
            % Syntax:
            %   slamSystem.setFixOlderPlatformVertices(unfixedTimeWindow)
            %
            % Description:
            %   Landmarks which are older than
            %   currentTime-unfixedTimeWindow have their status changed to
            %   fixed. This means they are not optimised. This saves both
            %   storage and computational costs.
            %
            % Inputs:
            %   unfixedTimeWindow - (double)
            %       The length of time a platform vertex remains unfixed
            obj.fixOlderPlatformVertices = true;
            obj.unfixedTimeWindow = unfixedTimeWindow;
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
            
            % Allocate the landmark map
            obj.landmarkIDStateVectorMap = containers.Map('KeyType', 'int64', 'ValueType', 'any');

            % Create the factor graph and set up the optimizer
            obj.graph = g2o.core.SparseOptimizer();
            %algorithm = g2o.core.GaussNewtonOptimizationAlgorithm();
            algorithm = g2o.core.LevenbergMarquardtOptimizationAlgorithm();
            obj.graph.setAlgorithm(algorithm);

            % Preallocate for convenience
            obj.platformVertices = cell(1, obj.INITIAL_CACHE_SIZE);
            
            % No vehicle vertices initally set
            obj.platformVertexId = 0;
            
            % The set of prediction edges, initially empty
            obj.processModelEdges = cell(1, obj.INITIAL_CACHE_SIZE);
            obj.numProcessModelEdges = 0;

            obj.currentPlatformVertex = {};
        end

        function stop(obj)
            % STOP Stop the SLAM system
            %
            % Syntax:
            %   slamSystem.stop()
            %
            % Description:
            %   Stop the SLAM system. Make sure all vertices are not fixed
            %   and run one final optimization pass.

            % Run the optimizer
            obj.optimize(20);


            % If we are fixing past vehicle states (Q3) then handle
            % unfixing for the final optimization pass
            if (obj.fixOlderPlatformVertices == true)

                % Get all the vertices apart from the first and make sure
                % they are not fixed. The first has to be fixed to embed
                % initial conditions.
                vertices = obj.graph.vertices();
                for v = 2 : numel(vertices)
                    vertices{v}.setFixed(false);
                end

                % Run the optimizer
                obj.optimize(50);
            end

            % Now clean up anything that needs cleaning up in the base
            % class.
            stop@cw1.drivebot.SLAMSystemBase(obj);
        end

        function chi2 = optimize(obj, maximumNumberOfOptimizationSteps)
            % OPTIMIZE Optimize the SLAM system
            %
            % Syntax:
            %   chi2 = slamSystem.optimize();
            %   slamSystem.optimize(maximumNumberOfOptimizationSteps);
            %
            % Description:
            %   Run the optimizer on g2o. Various performance metrics are
            %   collected as well. The optimizer is written so that you can
            %   call the method multiple times directly.
            %
            % Inputs:
            %   maximumNumberOfOptimizationSteps - (int)
            %       The maximum number of optimization steps to perform
            %       before giving up [Default: 10]
            %
            % Outputs:
            %   chi2 - (double)
            %       This is the sum of all the terms e'*Omega*e in the
            %       graph given the current converged estimate.

            % Record the time the optimization was called
            obj.performanceData.add('g2o.op.time', obj.currentTime);

            tStart = tic;
            obj.graph.initializeOptimization(obj.validateGraphOnInitialization);

            % Record initialization time
            obj.performanceData.add('g2o.op.io_dt', toc(tStart));

            % Set default number of optimization steps if not specified
            if (nargin == 1)
                maximumNumberOfOptimizationSteps = 10;
            end

            % Run the optimizer and record the time
            tStart = tic;                
            numIterations = obj.graph.optimize(maximumNumberOfOptimizationSteps);
            obj.performanceData.add('g2o.op.op_dt', toc(tStart));

            % Store the number of iterations
            obj.performanceData.add('g2o.op.it', numIterations);

            % Store the chi2 value
            chi2 = obj.graph.chi2();
            obj.performanceData.add('g2o.op.chi2', chi2);


            % Store vertex and edge sizes
            obj.performanceData.add('g2o.gr.n_ver', obj.graph.numVertices());         
            obj.performanceData.add('g2o.gr.n_edg', obj.graph.numEdges());         

            % Storage size for the optimization
            [b, H] = obj.graph.computeHessian();
            obj.performanceData.add('g2o.op.dim', numel(b));
            obj.performanceData.add('g2o.op.str', nnz(H) + numel(b));
        end

        function graph = optimizer(obj)
            % OPTIMIZER Return the underlying optimizer
            %
            % Syntax:
            %   graph = slamSystem.optimizer();
            %
            % Description:
            %   Every instance of this type of SLAM system has a single g2o
            %   graph. This method returns the instance to it.
            %
            % Outputs:
            %   graph - (g2o.core.SparseOptimizer)
            %       The graph object.

            graph = obj.graph;
        end

        function setValidateGraph(obj, validateGraphOnInitialization)
            % SETVALIDATEGRAPH Set the validate graph flag
            %
            % Syntax:
            %   slamSystem.setValidateGraph(validateGraphOnInitialization);
            %
            % Description:
            %   A properly configured graph obeys several conditions. These
            %   include all the vertex slots on all edges have been set,
            %   all the information matrices are positive semidefinite, and
            %   there are no straggling NaNs caused by values not being
            %   set. Furthermore, all vectors and matrices need to have the
            %   correct dimensions.
            %   
            %   If this flag is set to true, every edge and vertex is
            %   checked to make sure it is valid. However, this can greatly
            %   slow down the operation of the graph and should only be
            %   used for detailed debugging.
            %
            % Inputs:
            %   validateGraphOnInitialization - (bool)
            %       Flag to enable detailed checks.
            %
            % Outputs:
            %   chi2 - (double)
            %       This is the sum of all the terms e'*Omega*e in the
            %       graph given the current converged estimate.
            obj.validateGraphOnInitialization = validateGraphOnInitialization;
        end
        
        function validateGraphOnInitialization = validateGraph(obj)
            % VALIDATEGRAPH Return the validate graph flag
            %
            % Syntax:
            %   validateGraphOnInitialization = slamSystem.validateGraph();
            %
            % Description:
            %   Return the flag which determines if the graph should be
            %   validated when it is initialized.
            %
            %
            % Outputs:
            %   validateGraphOnInitialization - (bool)
            %       The value of the validate graph flag
            %
            % See also:
            %    SETVALIDATEGRAPH
            validateGraphOnInitialization = obj.validateGraphOnInitialization;
        end

        function [x,P] = platformEstimate(obj)
            % PLATFORMESTIMATE Return the mean and covariance of the platform estimate.
            %
            % Syntax:
            %   [x,P] = slamSystem.platformEstimate();
            %
            % Description:
            %   Run the optimizer and return the mean (platform vertex
            %   estimate) and covariance (approximate value computed from
            %   the Laplace approximation) of the platform state estimate.            %
            %
            % Outputs:
            %   x - (double vector)
            %       The value of the state
            %   P - (double PSD square matrix)
            %       The approximate covariance of the platform estimate.
            %
            % See also:
            %   LANDMARKESTIMATES

            obj.optimize(20);
            [x, P] = obj.graph.computeMarginals(obj.currentPlatformVertex);
            x=full(x{1});
            P=full(P{1});
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
            %   Note: this method does not call optimize. If you want to
            %   get the updated estimate, you must call optimize or
            %   retrieve the platform state estimate first.
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
            %
            % See also:
            %   PLATFORMESTIMATE

            landmarkVertices = values(obj.landmarkIDStateVectorMap);
            
            numberOfLandmarks = length(landmarkVertices);
            
            landmarkIds = NaN(1, numberOfLandmarks);
            m = NaN(obj.NL, numberOfLandmarks);
            Pmm = NaN(obj.NL, obj.NL, numberOfLandmarks);
            
            [x, P] = obj.graph.computeMarginals(landmarkVertices);

            % Now reformat. Grrr...
            for l = 1 : numberOfLandmarks
                m(:, l) = x{l};
                Pmm(:, :, l) = P{l};
            end            
        end

        
        function [T, X, PX] = platformEstimateHistory(obj)
            % PLATFORMESTIMATEHISTORY Deprecated method to return platform estimates over time.
            T = [];%obj.timeStore;
            X = [];%obj.xStore;
            PX = [];%obj.PStore;
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

            % Create the next vehicle vertex and add it to the graph; note
            % that the time is equal to the current time + dT to ensure
            % that it refers to "now"
            obj.currentPlatformVertex = ...
		cw1.drivebot.graph.PlatformStateVertex(obj.currentTime + dT);
            obj.graph.addVertex(obj.currentPlatformVertex);
            
            % Now create the process model edge and add it to the graph. We
            % can only do this if we have an odometry value
            processModelEdge = cw1.drivebot.graph.PlatformPredictionEdge(dT);
            processModelEdge.setVertex(1, obj.platformVertices{obj.platformVertexId});
            processModelEdge.setVertex(2, obj.currentPlatformVertex);
            processModelEdge.setMeasurement(obj.u);
            processModelEdge.setInformation(inv(obj.covU));
            processModelEdge.initialEstimate();
            obj.graph.addEdge(processModelEdge);
            
            obj.numProcessModelEdges = obj.numProcessModelEdges + 1;
            obj.processModelEdges{obj.numProcessModelEdges} = processModelEdge;
            
            % Bump the indices
            obj.platformVertexId = obj.platformVertexId + 1;
            obj.platformVertices{obj.platformVertexId} = obj.currentPlatformVertex;

            tStop = toc(tStart);

            % Record the time the optimization was called
            obj.performanceData.add('g2o.pr.time', obj.currentTime);
            obj.performanceData.add('g2o.pr.dt', tStop);


            % For Q3b, fix older platform vertices
            if (obj.fixOlderPlatformVertices == true)

                % Get the vertices
                vertices = obj.graph.vertices();

                % Figure out the freeze time - vertices from times older
                % than this are fixed
                fixTime = obj.currentTime - obj.unfixedTimeWindow;

                % Iterate over all platform state vertices and freeze the
                % ones which are too old
                for v = 1 : numel(vertices)
                    vertex = vertices{v};
                    if (contains(class(vertex), 'PlatformStateVertex') == true)
                        if (vertex.time() < fixTime)
                            vertex.setFixed(true);
                        end
                    end
                end
            end

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
            
            % Create the first vertex, set its estimate to the initial
            % value and add it to the graph.
            obj.currentPlatformVertex = cw1.drivebot.graph.PlatformStateVertex(obj.currentTime);
            obj.currentPlatformVertex.setEstimate(event.data);
            obj.graph.addVertex(obj.currentPlatformVertex);
            
            % Set the book keeping for this initial vertex.
            obj.platformVertexId = 1;
            obj.platformVertices{obj.platformVertexId} = obj.currentPlatformVertex;
            
            % If the covariance is 0, the vertex is known perfectly and so
            % we set it as fixed. If the covariance is non-zero, add a
            % unary initial prior condition edge instead. This adds a soft
            % constraint on where the state can be.
            if (det(event.covariance) < 1e-6)
                obj.currentPlatformVertex.setFixed(true);
            else
                initialPriorEdge = drivebot.graph.InitialPriorEdge();
                initialPriorEdge.setMeasurement(event.data);
                initialPriorEdge.setInformation(inv(event.covariance));
                initialPriorEdge.setVertex(obj.currentPlatformVertex);
                obj.graph.addEdge(initialPriorEdge);
            end

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

            obj.performanceData.add('g2o.slm.time', obj.currentTime);
            obj.performanceData.add('g2o.slm.nl', numel(event.info));

            tStart = tic;

            % Iterate over all the landmark measurements
            for l = 1 : numel(event.info)
                
                % Get the landmark vertex associated with this measurement.
                % If necessary, a new landmark vertex is created and added
                % to the graph.
                [landmarkVertex, newVertexCreated] = obj.createOrGetLandmark(event.info(l));
                z = event.data(:, l);
                
                % Add the measurement edge
                landmarkRangeBearingEdge = cw1.drivebot.graph.LandmarkRangeBearingEdge();
                landmarkRangeBearingEdge.setVertex(1, obj.currentPlatformVertex);
                landmarkRangeBearingEdge.setVertex(2, landmarkVertex);
                landmarkRangeBearingEdge.setMeasurement(z);
                landmarkRangeBearingEdge.setInformation(inv(event.covariance));
                
                if (newVertexCreated == true)
                    landmarkRangeBearingEdge.initialEstimate();
                end
                
                obj.graph.addEdge(landmarkRangeBearingEdge);

                % The code below supports Q3a. The number of edges a
                % landmark vertex particpates in tells you how many times
                % the landmark is being observed. If this exceeds a given
                % threshold we delete (or prune) excess edges. This is done
                % randomly here. The first measurement is always kept to ensure
                % that the map doesn't wander away over time. We also keep
                % the most recent measurement since that has just arrived.

                if (obj.limitObservationsPerLandmark == true)
                    if (landmarkVertex.numberOfEdges() > obj.maxObservationsPerLandmark)
                        edges = landmarkVertex.edges();
                        edgeToRemove = 1 + randi(obj.maxObservationsPerLandmark - 2);
                        obj.graph.removeEdge(edges{edgeToRemove});
                    end
                end
            end

            obj.performanceData.add('g2o.slm.dt', toc(tStart));

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

            % Create a GPS measurement edge
            gpsMeasurementEdge = cw1.drivebot.graph.GPSMeasurementEdge();
            gpsMeasurementEdge.setVertex(1, obj.currentPlatformVertex);
            gpsMeasurementEdge.setMeasurement(event.data);
            gpsMeasurementEdge.setInformation(inv(event.covariance));
            obj.graph.addEdge(gpsMeasurementEdge);

            success = true;
        end

        function success = handleCompassObservationEvent(obj, event)
            % HANDLECOMPASSOBSERVATIONEVENT Handle the compass observation event.
            %
            % Syntax:
            %   slamSystem.handleCompassObservationEvent(event);
            %
            % Description:
            %   Handle the compass observation.
            %
            % Inputs:
            %   event - (ebe.core.Event)
            %       The compass observation event.
            %
            % See Also:
            %   SIMULATOR

            % Add the edge
            compassMeasurementEdge = ...
                cw1.drivebot.graph.CompassMeasurementEdge(obj.configuration.compassAngularOffset);
            compassMeasurementEdge.setVertex(1, obj.currentPlatformVertex);
            compassMeasurementEdge.setMeasurement(event.data);
            compassMeasurementEdge.setInformation(inv(event.covariance));
            obj.graph.addEdge(compassMeasurementEdge);

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

            % Add the edge

            for s = 1 : numel(event.info)
                sensor = obj.map.sensors.bearing.sensors(event.info(s));
                bearingMeasurementEdge = ...
                    cw1.drivebot.graph.BearingMeasurementEdge(sensor.position, sensor.orientation);
                bearingMeasurementEdge.setVertex(1, obj.currentPlatformVertex);
                bearingMeasurementEdge.setMeasurement(event.data);
                bearingMeasurementEdge.setInformation(inv(event.covariance));
                obj.graph.addEdge(bearingMeasurementEdge);
            end

            success = false;
        end

        % This method returns a landmark associated with landmarkId. If a
        % landmark exists already, it is returned. If it does not exist, a
        % vertex is created and is added to the graph.
        function [landmarkVertex, newVertexCreated] = createOrGetLandmark(obj, landmarkId)
            
            % If the landmark exists already, return it
            if (isKey(obj.landmarkIDStateVectorMap, landmarkId) == true)
                landmarkVertex = obj.landmarkIDStateVectorMap(landmarkId);
                newVertexCreated = false;
                return
            end
            
            fprintf('Creating landmark %d\n', landmarkId);
            
            % Create the new landmark add it to the graph
            landmarkVertex = cw1.drivebot.graph.LandmarkStateVertex(landmarkId);
            obj.landmarkIDStateVectorMap(landmarkId) = landmarkVertex;
            
            obj.graph.addVertex(landmarkVertex);
            
            newVertexCreated = true;
        end

        function storeStepResults(obj)
            % % Store the estimate for the future
            % obj.timeStore(:, obj.stepNumber + 1) = obj.currentTime;
            % obj.xStore(:, obj.stepNumber + 1) = obj.x(1:cw1.drivebot.SystemModel.NP);
            % obj.PStore(:, obj.stepNumber + 1) = diag(obj.P(1:cw1.drivebot.SystemModel.NP, ...
            %     1:cw1.drivebot.SystemModel.NP));
        end
    end
end
