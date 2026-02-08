% Run q2b

import ebe.core.*;
import ebe.graphics.*;
import cw1.*;

% Find, load and parse the configuration file
config = ebe.utils.readJSONFile('config/q3_a.json');

% Create the mainloop object, which manages everything
mainLoop = ebe.MainLoop(config);

% Create the simulator and register it
simulator = drivebot.Simulator(config);
mainLoop.setEventGenerator(simulator);

% Create the G2O SLAM system and register it
g2oSLAMSystem = drivebot.G2OSLAMSystem(config);
mainLoop.addEstimator(g2oSLAMSystem);

g2oPrunedSLAMSystem = drivebot.G2OSLAMSystem(config);
g2oPrunedSLAMSystem.setMaxObservationsPerLandmark(4);
g2oPrunedSLAMSystem.setName('g2op-slam');
mainLoop.addEstimator(g2oPrunedSLAMSystem);

% Create the SLAM system and register it
ekfSLAMSystem = drivebot.EKFSLAMSystem(config);
mainLoop.addEstimator(ekfSLAMSystem);

% Collect performance data every 100 timesteps
mainLoop.setAccumulateResultsUpdatePeriod(25);

% Create the store for estimates
resultsAccumulator = ebe.slam.XPPlatformAccumulator();
mainLoop.addResultsAccumulator(resultsAccumulator);

% Set up the figure in which we draw everything
fig = FigureManager.getFigure("Q2c");
clf
hold on
axis([-5 55 -5 55])
axis equal

% Set up the views which show the output of the simulator
simulatorViewer = ebe.graphics.ViewManager(config);
simulatorView = drivebot.SimulatorView(config, simulator);
%simulatorView.setCentreAxesOnTruth(true);
simulatorViewer.addView(simulatorView);
simulatorViewer.addView(drivebot.SLAMSystemView(config, ekfSLAMSystem));
simulatorViewer.addView(drivebot.SLAMSystemView(config, g2oSLAMSystem));
simulatorViewer.addView(drivebot.SLAMSystemView(config, g2oPrunedSLAMSystem));

% Register the viewer with the mainloop
mainLoop.addViewer(simulatorViewer);

mainLoop.setGraphicsUpdatePeriod(25);

% Run the main loop until it terminates
mainLoop.run();

% Plot out state information

stateLabels = {'$x$','$y$', '$\theta$'};
axisLabels = {'Position', 'Position', 'Angle (rad)'};

TEstimator = resultsAccumulator.timeStore;
XTrueHistory = resultsAccumulator.xTrueStore;

for e = 1 : numel(resultsAccumulator.xEstStore)
    ebe.graphics.FigureManager.getFigure(sprintf('Results for %d', e));
    clf
    hold on

    PX = resultsAccumulator.PEstStore{e};
    X = resultsAccumulator.xEstStore{e};

    for f = 1 : 3
    
        subplot(3,1,f)
        sigmaBound = 2 * sqrt(PX(f, :));
        plot(TEstimator, -sigmaBound, 'r--', 'LineWidth', 2)
        hold on
        plot(TEstimator, sigmaBound, 'r--', 'LineWidth', 2)
        stateError = X(f, :) - XTrueHistory(f, :);
        if (f == 3)
            stateError = atan2(sin(stateError), cos(stateError));
        end
        plot(TEstimator, stateError, 'LineWidth', 2);
    
        % Work out the axes
        maxError = max(abs(stateError));
        bound = 1.1 * max(maxError, max(sigmaBound));
        axis([TEstimator(1) TEstimator(end) -bound bound])
        
        xlabel('Time (s)')
        ylabel(axisLabels{f}, 'Interpreter','latex')
        title(stateLabels{f}, 'Interpreter','latex')
    end
end

% Timing results
cw1.q3_a_show_results
