function sData = sample_connection_from_stable(stableChart, globalUnstableCoords, globalStableCoords, globalTime)
%SAMPLE_CONNECTION_FROM_STABLE - Given a candidate intersection for two charts, follow it backward through the chart lineage
% (forward in time) to the local stable manifold. Evaluation is returned in separate strands for each RegType.
%
%   Syntax:
%       output = SAMPLE_CONNECTION_FROM_STABLE(input)
%
%   Inputs:
%       input1 - Description
%       input2 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 23-Jan-2023;

%% Traverse the connection through the stable manifold
stableTime = globalTime(globalTime >= globalUnstableCoords(2)) - globalUnstableCoords(2); % global time points where the orbit lies in the stable atlas
stableTime = stableTime + globalStableCoords(2);
stableTime(end) = 0; % zero out final stable time to deal with rounding errors
sData.Orbit = {};
sData.Tau = {};
sData.RegVector = [];

% loop through stable lineage (backwards) starting from the intersection and moving forward in time until the local stable manifold is reached
chartStack =  stableChart.lineage(); % this is not flipped so we can move fowards from the intersection point through the stable lineage
regStack = [chartStack.RegType];
while ~isempty(regStack)
    nextStrandIdx = find_first_diff(regStack);
    strandIdx = 1:nextStrandIdx - 1;
    [strandOrbit, strandTau] = sample_strand(chartStack(strandIdx), stableTime, globalStableCoords(1));
    
    sData.RegVector(end+1) = chartStack(1).RegType; % append the regType of this orbit segment.
    sData.Orbit{end+1} = strandOrbit(:, 1:4).';  % end the previous orbit segment and append it to the orbit data array
    sData.Tau{end+1} = strandTau; % append the vector of shooting times for this segment
    
    chartStack = chartStack(nextStrandIdx:end);  % clip the charts for this strand off the stack
    regStack = regStack(nextStrandIdx:end);  % clip the vector of RegTypes to reflect what remains on the stack
end

% the initial point of the first strand is added when the strand is initialized. But this overshoots the point of
% intersection where the connection was found. So we throw it out.
sData.Orbit{1} = sData.Orbit{1}(:, 2:end);
sData.Tau{1} = sData.Tau{1}(:, [1, 3:end]);  % remove the connection time between the endpoint which was removed and the (now) first node
end % end sample_connection_from_stable




function [strandOrbit, strandTau] = sample_strand(chartVector, timeNodes,  globalSpatialCoordinate)
% SAMPLE_STRAND - sample a connecting orbit along a single strand of charts corresponding to a constant RegType

% Check if terminal point of first chart is already in evaluation nodes. If so, initialize empty arrays. If not, evaluate it
% to initialize arrays.
if ismember(chartVector(1).TimeSpan(2), timeNodes)  % terminal point of first chart is already in evaluation set
    strandTau = [];
    strandOrbit = [];
else  % evaluate the terminal point of first chart and initialize strand
    strandTau = 0;  % initialize a new vector of time nodes for this strand. Assume 0 for the first node which is the hotswap evaluation point.
    strandOrbit = chartVector(1).eval([globalSpatialCoordinate, 1], 'globalSpace', true); % initialize the array of evaluations at the nodes and evaluate the terminal point on the first chart
    strandOrbit = [strandOrbit{:}];  % unpack cell array into matrix
end
terminalOffset = 0;  % This keeps track of the global RegTime gap between the last evaluation of the previous chart in the strand and the terminal global regTime of the previous chart.
% In other words, it keeps track of how much time the current shooting segment spent in the previous chart so it can be deducted from the next chart in the strand.

% loop through the vector of charts for this strand and evaluate
for chart = chartVector
    [orbit, tau, terminalOffset] = sample_chart(chart, timeNodes, terminalOffset, globalSpatialCoordinate);
    strandOrbit = cat(1, strandOrbit, orbit); % append new orbit data to the top of the data
    strandTau = cat(2, strandTau, tau); % append RegTau vector for this chart
end

% Check if initial point of the last chart was already evaluated. If so, do nothing. If not, evaluate it and append it to the strand.
if ~ismember(chart.TimeSpan(1), timeNodes)  % initial point of last chart is not in evaluation set
    strandTau = cat(2, strandTau, terminalOffset);  % add terminal offset as the final Tau in RegTime
    terminalPt = chart.eval([globalSpatialCoordinate, 0], 'globalSpace', true);
    strandOrbit = cat(1, strandOrbit, [terminalPt{:}]);  % add initial evaluation to strand orbit
end
end


function [evals, tau, terminalOffset] = sample_chart(chart, fullNodes, terminalOffset, globalSpatialCoordinate)
% SAMPLE_CHART - sample a connecting orbit at given time nodes for a single Chart.

chartNodes = fullNodes(min(chart.TimeSpan) <= fullNodes & fullNodes <= max(chart.TimeSpan));  % identify all global F0 time nodes (if any) which lie in this chart's timespan
nChartNodes = length(chartNodes);  % count how many evaluation nodes lie in this chart

if nChartNodes > 0 % this chart has at least one evaluation node
    chartF0Time = chartNodes - chart.TimeSpan(1); % get global F0 time nodes for evaluations in this Chart
    chartRegTime = arrayfun(@(F0Tau)inverseregtime(chart, F0Tau), chartF0Time); % map time intervals for each shooting segment to regularized time as needed.
    chartLocalTime = chartRegTime./chart.LocalTime;  % map global RegTime to Local RegTime
    evalData = [globalSpatialCoordinate(1) * ones(nChartNodes, 1), chartLocalTime.']; % evaluation nodes are specified in local regTime coordinates and global space coordinates.
    cellEvals = chart.eval(evalData, 'globalSpace', true); % evaluation of this chart at the chartNodes.
    evals = [cellEvals{:}];
    tau = [terminalOffset + abs(chart.LocalTime - chartRegTime(1)), diff(chartRegTime)]; % append RegTime vector for evaluations
    terminalOffset = abs(chartRegTime(end)); % carry over the RegTime gap between the last chart evaluation and the chart's terminal RegTime
    
else  % no evaluations in this chart. Add its entire time span to the current offset
    evals = [];
    tau = [];
    terminalOffset = terminalOffset + abs(chart.LocalTime); % increment the running offset with the full Global RegTime interval for this chart
end
end

function idx = find_first_diff(arr)
% Find the least index where the specified vector changes value
idx = find(arr ~= arr(1), 1, 'first');
if isempty(idx)
    idx = length(arr) + 1;
end
end