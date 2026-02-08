% Simple static object example. This shows how to build a simple graph to
% estimate the state of a single object from a bunch of measurements of
% its position. This is a very simple version of using a GPS system to
% estimate the position of a static object.

import ebe.graphics.FigureManager;
import g2o.core.*;
import l4s.one_d_static.*;

% Some parameters
numberOfMeasurements = 100;
sigmaR(1) = 1;
sigmaR(2) = 10;

% The information (inverse covariance) for the measurement edge
omegaR = 1 ./ sigmaR;

% Ground truth location
trueX = 10;

% Odd even flag. This causes us to "switch" between the two sensors. The
% first line means we always use the good sensor, the second case
% alternates between the two
%oddEven = ones(1, numberOfMeasurements);
oddEven = mod(0:numberOfMeasurements-1, 2) + 1;

% Sample the noises for the different observations
z = trueX + sqrt(sigmaR(oddEven)) .* randn(numberOfMeasurements, 1);

% Create the graph
graph = SparseOptimizer();
algorithm = GaussNewtonOptimizationAlgorithm();
graph.setAlgorithm(algorithm);

% Create the vertex. This contains the state we want to estimate.
v = ParticleStateVertex();

% Set the initial estimate. We have to set some initial value. This ideally
% shouldn't be too far from the final soluton. Here we use the first
% measurement.
v.setEstimate(z(1));

% Added the vertex to the graph. The graph now knows that we have some
% states we want to estimate.
graph.addVertex(v);

% Store for the covariance and the estimate
xStore = NaN(1, numberOfMeasurements);
PStore = NaN(1, numberOfMeasurements);

% Create the edges. Each edge corresponds to a position measurement.

for k = 1 : numberOfMeasurements

    % Create the measurement edge
    e = ParticleMeasurementEdge();

    % Link it so that it connects to the vertex we want to estimate
    e.setVertex(1, v);

    % Set the measurement value and the measurement covariance
    e.setMeasurement(z(k));
    e.setInformation(omegaR(oddEven(k)));

    % Add the edge to the graph; the graph now knows we have these edges
    % which need to be added
    graph.addEdge(e);

    % Graph construction complete for this iteration

    % Initialise the optimization. If you ever change the structure of the
    % graph (e.g., add remove edges and vertices) you must call this before
    % calling optimize. If you don't change the graph structure, you can call
    % optimze multiple times without calling initializeOptimization().
    graph.initializeOptimization();

    % Run the optimizer. By default it does at most 10 iterations, but you can
    % pass the number in optionally.
    graph.optimize();

    % Get the estimate and covariance. In graph-land, getting the covariance
    % matrix is described as "computing the marginals". Note the vector X is
    % the "big state vector" which consists of the states from all the vertices
    % concatenated on one another. Similarly, PX is a sparse covariance matrix
    % which contains the covariances of the vertices.
    [X, PX] = graph.computeMarginals();
    
    % Store the estimate and covariance. In this case, because they are
    % scalar, we simply copy them across.
    xStore(k) = X{1};
    PStore(k) = PX{1};
end

% Plot out state information
ebe.graphics.FigureManager.getFigure('Lab 04 Activity 01 State Error');
clf
hold on
subplot(3,1,1)
plot(xStore)
xlabel('Measurement number')
ylabel('Estimate')
subplot(3,1,2)
plot(sqrt(PStore))
xlabel('Measurement number')
ylabel('Standard deviation')
subplot(3,1,3)
hold on
plot(xStore-trueX)
plot(2*sqrt(PStore),'r--')
plot(-2*sqrt(PStore),'r--')
xlabel('Measurement number')
ylabel('Error')