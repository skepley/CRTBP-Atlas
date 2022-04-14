function connection = connection2BVP(connectionData, nNode, stableLocalMap, unstableLocalMap)
%CONNECTION2BVP - Process a connecting orbit from the mining algorithm into an approximation for the BVP solver.
%   CONNECTION2BVP() - A more detailed description of the function
%
%   Syntax:
%       output = CONNECTION2BVP(input1, input2)
%       [output1, output2] = CONNECTION2BVP(input1, input2, input3)
%
%   Inputs:
%       connectionData - Cell array of the form: {stableChart, unstableChart, connectionCoordinates}
%           Both Charts are RegCRTBPChart objects which parameterized invariant manifolds and the connection 
%           Coordinates are a vector in R^3 of the form: (s1, t1, s2) such that stable(s1, t1) = unstable(s2, 0)
%       nNode - Number of uniformly spaced time points to evaluate the orbit. These evaluations are the endpoints of the shooting segments
%       stableLocalMap - A function handle which returns the preimage of the local stable manifold parameterization at the endpoint where 
%           it meets the connecting orbit.
%       unstableLocalMap - A function handle which returns the preimage of the local unstable manifold parameterization at the endpoint where 
%           it meets the connecting orbit.
%
%   Outputs:
%       connection - A struct with the data necessary to form an inital BVP solution guess. 
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 20-Aug-2020; Last revision: 20-Aug-2020



%% ================================================== PARSE CONNECTION PROPERTIES ==================================================
warning('OFF', 'Chart:eval');  % turn off Chart.eval warning messages
warning('OFF', 'Chart:local2global');  % turn off Chart.eval warning messages
stableChart = connectionData{1};
unstableChart = connectionData{2};
localIntersection = connectionData{3};  % x = (s1, t1, s2) and S(s1, t1) = U(s2, t0) in 3 coordinates
localStableSpace = localIntersection(1); % local spatial coordinate for the stable chart at the intersection
localStableRegTau = localIntersection(2); % local time for the stable chart at the intersection
localStableF0Tau = taylorregtime(stableChart, localStableRegTau);
localUnstableSpace = localIntersection(3); % local spatial coordinate for the unstable chart at the intersection
localUnstableF0Tau = 0; % regularize local time with t0 = 0.

globalStableCoords = stableChart.local2global([localStableSpace, localStableF0Tau]);  % global intersection coordinates (Ss, Ts)
globalUnstableCoords = unstableChart.local2global([localUnstableSpace, localUnstableF0Tau]);  % global intersection coordinates (Us, Ut)
globalConnectionTime = globalUnstableCoords(2) - globalStableCoords(2);  % time of flight in F0 time
globalTime = linspace(0, globalConnectionTime, nNode); % a uniform grid of global F0 time points along the entire connecting orbit

%% ================================================== GENERATE THE CONNECTING ORBIT DATA ==================================================
connection = struct; % initialize the connection struct
connection.ConnectionTime = globalConnectionTime;  % time of flight from one local manifold to the other in F0 time units
connection.UnstableTime = globalUnstableCoords(2);  % time of flight captured in the unstable manifold
% map the local manifold data
connection.LocalUnstable = unstableLocalMap(globalUnstableCoords(1));
connection.LocalStable = stableLocalMap(globalStableCoords(1));
connection.RegVector = []; % initialize vector of segment regTypes
connection.Orbit = {};  % initialize cell array of orbit segments distinguished by regType
connection.Tau = {};  % initialize cell array of shooting times along orbit segments distinguished by regType

% ================================================== GENERATE UNSTABLE ORBIT SEGMENTS ==================================================
unstableTime = globalTime(globalTime < globalUnstableCoords(2));  % global time points where the orbit lies in the unstable atlas
% lineage map
unstableLineage = unstableChart;  % initialize unstable lineage as a RegCRTBPChart type and and empty it
unstableLineage = unstableLineage(2:end);
thisChart = unstableChart;
while ~isempty(thisChart)
    unstableLineage(end+1) = thisChart;
    thisChart = thisChart.ParentHandle;
end
unstableLineage = flip(unstableLineage);  % this is flipped so we can move backwards through the unstable lineage

% iniitalize first unstable orbit segment
orbitSegment = [];  % initialize empty segment under the first chart regularization type
tauSegment = []; % initialize empty segment of timesteps. Each entry is the time of flight between points on the orbit segment
terminalRegTauOffset = 0;
segmentRegType = unstableLineage(1).RegType;    % RegType of initial orbit segment

% loop over the unstable charts and evaluate the orbit points and shooting times
for j = 1:length(unstableLineage)
    jChart = unstableLineage(j);
    chartGlobalTime = unstableTime(min(jChart.TimeSpan) <= unstableTime & unstableTime < max(jChart.TimeSpan));  % identify all F0 time points which lie in this chart
    
    if ~isequal(jChart.RegType, segmentRegType) % this is the first chart of a new orbit segment
        
        % append the terminal point of the last chart to this orbit segment
        tauSegment = cat(2, tauSegment, terminalRegTauOffset); % append regTau for the additional evaluation node
        evalData = [globalUnstableCoords(1), 1]; % evaluate at local time with global space
        jOrbitData = unstableLineage(j-1).eval(evalData, 'globalSpace', true); % evaluation after switching to local coordinates
        orbitSegment = cat(1, orbitSegment, [jOrbitData{:}]); % append new orbit data to the top of the data
        
        % terminate the orbit segment
        tauSegment(1) = 0;  % zero out first entry to deal with rounding errors
        connection.RegVector(end+1) = segmentRegType; % append the regType of this orbit segment. 
        connection.Orbit{end+1} = orbitSegment(:, 1:4).';  % end the previous orbit segment and append it to the orbit data array
        connection.Tau{end+1} = tauSegment; % append the vector of shooting times for this segment
        
        % start a new orbit segment for this regType
        segmentRegType = jChart.RegType;  % update segment regtype variable
        tauSegment = [];  % initialize a new vector of shooting times for this regType with an assumed 0 for the first node which is the hotswap point
        orbitSegment = []; % initialize a new orbit segment
        terminalRegTauOffset = 0;
        chartGlobalTime = sort([chartGlobalTime, jChart.TimeSpan(1)]);  % add an initial evaluation node to the first chart
    end
      
    if ~isempty(chartGlobalTime) % this chart has at least one evaluation point
        chartF0Tau = chartGlobalTime - jChart.TimeSpan(1); % get F0 global timesteps for this chart
        chartRegTau = arrayfun(@(F0Tau)inverseregtime(jChart, F0Tau), chartF0Tau); % map time intervals for each shooting segment to F1/F2 time as needed.
        nChartNode = length(chartGlobalTime);  % count how many evaluation nodes lie in this chart
        chartLocalTau = chartRegTau./jChart.LocalTime;  % map global RegTau to Local RegTau
        evalData = [globalUnstableCoords(1) * ones(nChartNode, 1), chartLocalTau.']; % evaluate at local time with global space
        jOrbitData = jChart.eval(evalData, 'globalSpace', true); % evaluation after switching to local coordinates
        orbitSegment = cat(1, orbitSegment, [jOrbitData{:}]); % append new orbit data to the top of the data
        tauSegment = cat(2, tauSegment, [terminalRegTauOffset + chartRegTau(1), diff(chartRegTau)]); % append RegTau vector for this chart
        terminalRegTauOffset = jChart.LocalTime - chartRegTau(end); % carry over the RegTime gap between the last chart evaluation and its terminal time
        
    else
        terminalRegTauOffset = terminalRegTauOffset + jChart.LocalTime; % increment the running offset with the full Global RegTime interval for this chart
    end
end
connectionOffset = globalUnstableCoords(2) - unstableTime(end);
terminalRegTauOffset = inverseregtime(unstableLineage(end-1), connectionOffset); % update Tau offset

% DEBUGGING THE STABLE/UNSTABLE ORBIT SEGMENTS
connection.Unstable = connection.Orbit;  % copy the current orbit segments so I can plot only half of the connection
connection.UnstableTau = connection.Tau;
connection.Unstable{end+1} = orbitSegment(:, 1:4).';  % add the unstable data from the current orbit segment
connection.UnstableTau{end+1} = tauSegment; 

%% ================================================== GET STABLE SEGMENTS ==================================================
stableTime = globalTime(globalTime >= globalUnstableCoords(2)) - globalUnstableCoords(2); % global time points where the orbit lies in the unstable atlas
stableTime = stableTime + globalStableCoords(2);
stableTime(end) = 0; % zero out final stable time to deal with rounding errors

% This switches the timestep order since the stable chart lineage is looped through backwards instead of backwards.
% Subtraciting the unstable connection time ensures that 0 will be an evaluation time i.e. the point on the local stable manifold is
% the final point on the orbit.
% lineage map
stableLineage = stableChart;  % initialize stable lineage as a RegCRTBPChart type and and empty it
stableLineage = stableLineage(2:end);
thisChart = stableChart;
while ~isempty(thisChart)
    stableLineage(end+1) = thisChart;
    thisChart = thisChart.ParentHandle;
end
% We do not flip the order of the stableLineage since we start evaluating at the intersection and work backward through the atlas.


% initialize stable part of orbit. The mining algorithm ensures that the stable/unstable intersection charts have the same
% RegType. So the first chart of the stable lineage will continue appending to the last orbit segment of the unstable lineage.
connectionOffset = stableLineage(1).TimeSpan(2) - globalStableCoords(2);  % adjust terminal chart time to begin at the connection time
terminalRegTauOffset = terminalRegTauOffset + inverseregtime(stableLineage(1), connectionOffset);

% loop over the stable charts and evaluate the orbit points and shooting times
for j = 1:length(stableLineage)
    jChart = stableLineage(j);
    chartGlobalTime = stableTime(min(jChart.TimeSpan) < stableTime & stableTime <= max(jChart.TimeSpan));  % identify all F0 time points which lie in this chart
    
    if ~isequal(jChart.RegType, segmentRegType) % this is the first chart of a new orbit segment
        
        % append the initial (since we are moving backward through the chart lineage) point of the previous chart to this orbit segment
        tauSegment = cat(2, tauSegment, terminalRegTauOffset); % append regTau for the additional evaluation node
        evalData = [globalStableCoords(1), 0]; % evaluate at local initial time with global space
        jOrbitData = stableLineage(j-1).eval(evalData, 'globalSpace', true); % evaluation after switching to local coordinates
        orbitSegment = cat(1, orbitSegment, [jOrbitData{:}]); % append new orbit data to the top of the data
        
        % terminate the orbit segment
        tauSegment(1) = 0;  % zero out first entry to deal with rounding errors
        connection.RegVector(end+1) = segmentRegType; % append the regType of this orbit segment. 
        connection.Orbit{end+1} = orbitSegment(:, 1:4).';  % end the previous orbit segment and append it to the orbit data array
        connection.Tau{end+1} = tauSegment; % append the vector of shooting times for this segment
        
        % start a new orbit segment for this regType
        segmentRegType = jChart.RegType;  % update segment regtype variable
        tauSegment = [];  % initialize a new vector of shooting times for this regType with an assumed 0 for the first node which is the hotswap point
        orbitSegment = []; % initialize a new orbit segment
        terminalRegTauOffset = 0;
        chartGlobalTime = sort([chartGlobalTime, jChart.TimeSpan(2)]);  % add a terminal evaluation node to the first chart
    end
    
    if ~isempty(chartGlobalTime) % this chart has at least one evaluation point
        chartF0Tau = chartGlobalTime - jChart.TimeSpan(1); % get F0 global timesteps for this chart
        chartRegTau = arrayfun(@(F0Tau)inverseregtime(jChart, F0Tau), chartF0Tau); % map F0 time intervals for each shooting segment to RegTime as needed.
        nChartNode = length(chartGlobalTime);  % count how many evaluation nodes lie in this chart
        chartLocalTau = chartRegTau./jChart.LocalTime;  % map global RegTau to Local RegTau
        evalData = [globalStableCoords(1) * ones(nChartNode, 1), chartLocalTau.']; % evaluate at local time with global space
        jOrbitData = jChart.eval(evalData, 'globalSpace', true); % evaluation after switching to local coordinates
        orbitSegment = cat(1, orbitSegment, [jOrbitData{:}]); % append new orbit data to the top of the data
        tauSegment = cat(2, tauSegment, [terminalRegTauOffset + abs(jChart.LocalTime - chartRegTau(1)), diff(chartRegTau)]); % append RegTau vector for this chart
        terminalRegTauOffset = abs(chartRegTau(end)); % carry over the RegTime gap between the last chart evaluation and its terminal time
        
    else
        terminalRegTauOffset = terminalRegTauOffset + abs(jChart.LocalTime); % increment the running offset with the full Global RegTime interval for this chart
    end
end
% write the final orbit segment
tauSegment(1) = 0;  % zero out first entry to deal with rounding errors
connection.RegVector(end+1) = segmentRegType; % append the regType of this orbit segment.
connection.Orbit{end+1} = orbitSegment(:, 1:4).';  % end the finaal orbit segment and append it to the orbit data array
connection.Tau{end+1} = tauSegment; % append the vector of shooting times for this segment
end % end connection2BVP

% Revision History:
%{

%}
