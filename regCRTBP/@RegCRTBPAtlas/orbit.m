function orbitData = orbit(obj, globalSpace, globalTime, varargin)
%ORBIT - Evaluate an orbit through the atlas
%
%   ORBIT() - A more detailed description of the function
%
%   Syntax:
%       orbitData = ORBIT(obj, s, T) returns an array of evaluations of a trajectory through the point x0(s) at the timesteps specified in T
%
%   Inputs:
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
%   Date: 23-Mar-2019; Last revision: 23-Mar-2019

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
    orbitData = cat(1, [thisOrbitData{:}], orbitData); % append new orbit data to the top of the data
    thisChart = thisChart.ParentHandle; % continue to next chart containing this orbit
end
end % end orbit

% Revision History:
%{

%}
