% This script shows how to set up the individual pieces of the code and
% start running them

import ebe.core.*;
import ebe.graphics.*;
import l1s.pointbot.*;

% Find, load and parse the configuration file
config = ebe.utils.readJSONFile('config/activity3.json');

% Set up the simulator and start it
simulator = Simulator(config);
simulator.start();

% Create the Kalman filter and start it
kf = KalmanFilter(config);
kf.start();

% Set up the figure in which we draw everything
fig = FigureManager.getFigure("Simulator Output");
clf
hold on
axis([-20 20 -20 20])
axis square

% Set up the view which show the output of the simulator
simulatorViewer = ebe.graphics.ViewManager(config);
simulatorViewer.addView(SimulatorView(config, simulator));

% Register a viewer to show the Kalman filter results
simulatorViewer.addView(ebe.graphics.MeanCovarianceView(config, kf, [1 3]));

simulatorViewer.start();

% This is the main loop - while the simulator says we should keep going, we
% step the simulator, extract the events which were generated, and
% visualize them.
while (simulator.keepRunning()  == true)
    simulator.step();
    events = simulator.events();
    kf.processEvents(events);
    simulatorViewer.visualize(events);
    drawnow
end
