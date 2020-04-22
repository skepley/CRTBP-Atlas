function varargout = patch(obj, gridArray, plotIdx, varargin)
%PATCH - Plot an atlas of Charts as a patch
%
%   PATCH() - A more detailed description of the function
%
%   Syntax:
%       output = PATCH(input1, input2)
%       [output1, output2] = PATCH(input1, input2, input3)
%
%   Inputs:
%       obj - An atlas of Chart objects parameterized on a subset of [-1,1] x [0, Tau]
%       gridArray - A cell array of global space and time evaluation nodes
%       idx - coordinate indices to plot
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
%   email: shane.kepley@rutgers.edu
%   Date: 19-Apr-2019; Last revision: 2-Nov-2019

% TODO:
% 1. Pass patch plotting options through as varargs

globalSpace = gridArray{1};
globalTime = gridArray{2};

atlasDimension = length(gridArray);
if ~isequal(atlasDimension, 2)
    error('not implemented')
end

% initialize patch data
face = [];
vertex = [];
colorData = [];

% get space/time boundary data
minSpace = min(globalSpace);
maxSpace = max(globalSpace);
minTime = min(globalTime);
maxTime = max(globalTime);

% main loop
for iChart = obj.Chart
    % get global evaluation node data for this chart
    
    % filter out global time coordinates
    t0 = min(iChart.TimeSpan);
    t1 = max(iChart.TimeSpan);
    iGlobalTimeIdx = (t0 <= globalTime) & (globalTime <= t1); % check data for evaluations which lie in this chart
    iTimeGridData = globalTime(iGlobalTimeIdx); % filter out valid global time evaluations
    
    % append boundary time nodes if this chart is in the interior of the time evaluation domain
    if t0 > minTime
        iTimeGridData = [t0, iTimeGridData];
    end
    if t1 < maxTime
        iTimeGridData = [iTimeGridData, t1];
    end
    
    % filter out global space coordinates
    s0 = iChart.SpatialSpan(1);
    s1 = iChart.SpatialSpan(2);
    iGlobalSpaceIdx = (s0 <= globalSpace) & (globalSpace <= s1); % check data for evaluations which lie in this chart
    iSpaceGridData = globalSpace(iGlobalSpaceIdx); % filter out valid global time evaluations
    
    % append boundary space nodes if this chart is in the interior of the spatial evaluation domain
    if s0 > minSpace
        iSpaceGridData = [s0, iSpaceGridData];
    end
    if s1 < maxSpace
        iSpaceGridData = [iSpaceGridData, s1];
    end
    
    
    % build a meshgrid from local evaluation grids
    [sMesh, tMesh] = meshgrid(iSpaceGridData, iTimeGridData);
    iEvalData = [reshape(sMesh, [], 1), reshape(tMesh, [], 1)];
    
    try
        iFace = delaunay(iEvalData(:,1), iEvalData(:,2)); % index triples for delaunay triangulation in space-time.
 
        % get global vertex data for this chart
        iPatchCell = iChart.eval(iEvalData, 'GlobalTime', true, 'GlobalSpace', true); % cell array of evaluations
        iPatchEval = cell2mat(iPatchCell(plotIdx)); % convert to data array consisting only of indices to plot
        iVertex = mid(iPatchEval); % convert to floats if necessary
        vertexCount = size(vertex, 1); % current vertex count to shift global vertex indices of this chart
        vertex = [vertex;iVertex]; % append new vertices
        
        % get global face data for this chart
        colorData = [colorData; iEvalData(:,2)]; % append global time data for patch coloring
        face = [face;iFace + vertexCount]; % start counting face labels from the previous end label
    catch % exception handling to deal with bad calls to delaunay
        
    end
end

% plot it
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData, 'FaceColor', 'interp', 'EdgeColor', 'none');
end % end patch

% Revision History:
%{
2-Nov-2019 - Fixed a bug which caused it to incorrectly attach charts to one another. Also changed the format so it takes
        an array of independent grids and computes a meshgrid for each individual Chart instead of computing a single
        meshgrid and indexing into it with each Chart. The pictures look much smoother this way.
%}
