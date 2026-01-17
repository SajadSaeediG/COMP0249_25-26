% This task demonstrates the simulator, using the MainLoop class, which
% wraps everything up.

import ebe.core.*;
import ebe.graphics.*;
import l1s.pointbot.*;

% Find, load and parse the configuration file
config = ebe.utils.readJSONFile('config/activity5.json');

% Create the mainloop object, which manages everything
mainLoop = ebe.MainLoop(config);

% Create the simulator and register it
simulator = Simulator(config);
mainLoop.setEventGenerator(simulator);

% Create the Kalman filter and register it
kf = KalmanFilter(config);
mainLoop.addEstimator(kf);

% Set up the figure in which we draw everything
fig = FigureManager.getFigure("Scenario Output");
clf
hold on
axis([-30 30 -30 30])
axis square

% Set up the views which show the output of the simulator
simulatorViewer = ebe.graphics.ViewManager(config);
simulatorViewer.addView(SimulatorView(config, simulator));

% Register a viewer to show the Kalman filter results
simulatorViewer.addView(ebe.graphics.MeanCovarianceView(config, kf, [1 3]));

% Register the viewer with the mainloop
mainLoop.addViewer(simulatorViewer);

% Run the main loop until it terminates
mainLoop.run();

% Now extract the estimate history from the estimation algorithm
[TEstimator, X, PX] = kf.estimateHistory();

[TSimulator, XTrueHistory] = simulator.history();

% Plot out the estimation errors and two standard deviation bounds
ebe.graphics.FigureManager.getFigure('Error Results');
clf

stateLabels = {'$x$','$\dot{x}$','$y$','$\dot{y}$'};

for f = 1 : 4
    subplot(4,1,f)
    sigmaBound = 2 * sqrt(PX(f, :));
    plot(TEstimator, -sigmaBound, 'r--', 'LineWidth', 2)
    hold on
    plot(TEstimator, sigmaBound, 'r--', 'LineWidth', 2)
    stateError = X(f, :) - XTrueHistory(f, 1:size(X,2));
    plot(TEstimator, stateError, 'b', 'LineWidth', 2);
    % Work out the axes
    maxError = max(abs(stateError));
    bound = 1.1 * max(maxError, max(sigmaBound));
    axis([TEstimator(1) TEstimator(end) -bound bound])

    xlabel('Time (s)')
    title(stateLabels{f}, 'Interpreter','latex')
    if (rem(f, 2) == 0)
        ylabel('Velocity $(ms^{-1})$', 'Interpreter','latex')
    else
        ylabel('Position $(ms)$', 'Interpreter','latex')
    end
    title('Results')
end