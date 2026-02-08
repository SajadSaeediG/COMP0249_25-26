% Process results

clear optimTimes chi2Values optimDurations

g2oPerfData = g2oSLAMSystem.getPerformanceData();
optimTimes{1} = g2oPerfData.get('g2o.op.time');
chi2Values{1} = g2oPerfData.get('g2o.op.chi2');
optimDurations{1} = g2oPerfData.get('g2o.op.op_dt');

% Graph stats
graph = g2oSLAMSystem.graph();

% Count the elements in the graph
hypergraphElements = {ebe.utils.PerformanceData()};
vertices = graph.vertices();
for v = 1 : length(vertices)
    hypergraphElements{1}.add(class(vertices{v}), 1);
end
edges = graph.edges();
for e = 1 : length(edges)
    hypergraphElements{1}.add(class(edges{e}), 1);
end

g2oPerfData = g2oPrunedSLAMSystem.getPerformanceData();
optimTimes{2} = g2oPerfData.get('g2o.op.time');
chi2Values{2} = g2oPerfData.get('g2o.op.chi2');
optimDurations{2} = g2oPerfData.get('g2o.op.op_dt');

graph = g2oPrunedSLAMSystem.graph();

% Count the elements in the pruned graph
hypergraphElements{2} = ebe.utils.PerformanceData();
vertices = graph.vertices();
for v = 1 : length(vertices)
    hypergraphElements{2}.add(class(vertices{v}), 1);
end
edges = graph.edges();
for e = 1 : length(edges)
    hypergraphElements{2}.add(class(edges{e}), 1);
end

fig = ebe.graphics.FigureManager.getFigure("Q3a Timing and chi2 results");
clf
subplot(2,1,1)
plot(chi2Values{1})
hold on
plot(chi2Values{2})
xlabel('Time step')
ylabel('chi2 values')
subplot(2,1,2)
plot(optimDurations{1}, 'g')
hold on
plot(optimDurations{2}, 'r')
xlabel('Time step')
ylabel('Time to optimize')
hold on
subplot(2,1,1)
