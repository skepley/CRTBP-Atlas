function orbitData = CRTBPorbit(obj, globalSpace, globalTime, varargin)
%CRTBPORBIT - Overloads the orbit method for evaluating an orbit in an atlas of RegCRTBPChart objects.
%
%   CRTBPORBIT() - This deals with the problem of the vector field dimensions changing depending on which coordinates are being used for a chart. The orbit
%       returned is always the first 4 coordinates with respect to F0 coordinates.
%
%   Syntax:
%       orbitData = ORBIT(obj, s, T) returns an array of evaluations of a trajectory through the point x0(s) at the timesteps specified in T
%
%   Inputs:
%       obj - An atlas whose Charts are RegCRTBPChart objects
%       globalSpace - A single double in the interval [-1,1]
%       globalTime - A vector of real floats at which to evaluate the trajectory
%
%   Outputs:
%       orbitData - An array of evaluations along the orbit
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Jul-2019; Last revision: 04-Jul-2019


%% find the terminal chart for this orbit
generationIdx = obj.LastGeneration; % start searching atlas at the last generation
while ~exist('terminalChart')
    thisGeneration = [obj.Chart([obj.Chart.Generation] == generationIdx)]; % get list of charts for this generation
    % check if orbit terminates in this generation
    j = 1;
    
    while j <= length(thisGeneration)
        jChart = thisGeneration(j);
        chkSpace = jChart.local2global(-1, 1) <= globalSpace && globalSpace <= jChart.local2global(1, 1); % check if this chart contains part of the orbit
        if chkSpace
            terminalChart = jChart; % found the terminal chart for this orbit
            j = length(thisGeneration) + 1; % break while loop
        else
            j = j+1;
        end
    end
    generationIdx = generationIdx-1; % check previous generation for terminal chart
end

%% loop through lineage of terminal chart and evaluate orbit
orbitData = []; % initialize orbit
thisChart = terminalChart;

while ~isempty(thisChart)
    globalTime = reshape(globalTime, [],1); % ensure time is a column vector
    thisGlobalData = [globalSpace*ones(size(globalTime)), globalTime]; % evaluation data in global coordinates
    thisOrbitData = thisChart.eval(thisGlobalData, 'globalTime', true, 'globalSpace', true); % evaluation after switching to local coordinates
    if ~isempty(thisOrbitData)
        % swap to F0 field
        switch thisChart.RegType
            case 0
                orbitSegment = [thisOrbitData{1:4}];
            case 1
                orbitSegment = CRTBP2reg([thisOrbitData{1:4}], thisChart.Parameter(1), -1);
            case 2
                orbitSegment = CRTBP2reg([thisOrbitData{1:4}], thisChart.Parameter(1), -2);
        end % switch
        orbitData = cat(1, orbitSegment(:, 1:4), orbitData); % append new orbit ****FROM THE STATE VARIABLES ONLY**** data to the top of the data
    end
    thisChart = thisChart.ParentHandle; % continue to next chart containing this orbit
end
end % end CRTBPorbit

% Revision History:
%{

%}
