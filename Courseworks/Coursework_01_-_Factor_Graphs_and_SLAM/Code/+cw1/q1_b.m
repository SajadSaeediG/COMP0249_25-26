% This script runs an intro scenarioinitializeop

import ebe.core.*;
import ebe.graphics.*;
import cw1.*;

% Find, load and parse the configuration file
config = ebe.utils.readJSONFile('config/q1_b.json');

% Create the mainloop object, which manages everything
mainLoop = ebe.MainLoop(config);

% Create the simulator and register it
simulator = drivebot.Simulator(config);
mainLoop.setEventGenerator(simulator);

% Create the SLAM system and register it
g2oSLAMSystem = drivebot.G2OSLAMSystem(config);
mainLoop.addEstimator(g2oSLAMSystem);

% Create the store for estimates
resultsAccumulator = ebe.slam.XPPlatformAccumulator();
mainLoop.addResultsAccumulator(resultsAccumulator);
mainLoop.setAccumulateResultsUpdatePeriod(50);

% Set up the figure in which we draw everything
fig = FigureManager.getFigure("Q1b");
clf
hold on
axis([-25 25 -25 25])
axis square

% Set up the views which show the output of the simulator
simulatorViewer = ebe.graphics.ViewManager(config);
simulatorView = drivebot.SimulatorView(config, simulator);
simulatorView.setCentreAxesOnTruth(true);
simulatorViewer.addView(simulatorView);
simulatorViewer.addView(drivebot.SLAMSystemView(config, g2oSLAMSystem));

% Register the viewer with the mainloop
mainLoop.addViewer(simulatorViewer);
mainLoop.setGraphicsUpdatePeriod(50);

% Run the main loop until it terminates
mainLoop.run();

% Plot out state information

stateLabels = {'$x$','$y$'};

TEstimator = resultsAccumulator.timeStore;
XTrueHistory = resultsAccumulator.xTrueStore;

for e = 1 : numel(resultsAccumulator.xEstStore)
    ebe.graphics.FigureManager.getFigure(sprintf('Results for %d', e));
    clf
    hold on

    PX = resultsAccumulator.PEstStore{e};
    X = resultsAccumulator.xEstStore{e};

    for f = 1 : 2
    
        subplot(2,1,f)
        sigmaBound = 2 * sqrt(PX(f, :));
        plot(TEstimator, -sigmaBound, 'r--', 'LineWidth', 2)
        hold on
        plot(TEstimator, sigmaBound, 'r--', 'LineWidth', 2)
        stateError = X(f, :) - XTrueHistory(f, :);
        plot(TEstimator, stateError, 'LineWidth', 2);
    
        % Work out the axes
        maxError = max(abs(stateError)) + 0.1;
        bound = 1.1 * max(maxError, max(sigmaBound));
        axis([TEstimator(1) TEstimator(end) -bound bound])
        
        xlabel('Time (s)')
        ylabel('Position $(ms)$', 'Interpreter','latex')
        title(stateLabels{f}, 'Interpreter','latex')
    end
end