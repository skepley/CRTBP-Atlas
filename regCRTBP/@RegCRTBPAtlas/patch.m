function varargout = patch(obj, evalNode, chartRegType, plotRegType, plotIdx, varargin)
%PATCH - Modified version of the Atlas patch method which includes hotswapping
%
%   PATCH() - A more detailed description of the function
%
%   Syntax:
%      PATCH(obj, {S,T}, R1, R2, plotIdx) plots the charts in the RegCRTBPAtlas (obj) which are computed with respect to 
%       R1 (0, 1, or 2) regtype using coordinates with respect to the R2 regtype at the rectangular space/time mesh with nodes S and T
%       respectively. plotIdx is a 2 or 3 vector which identifies which of the 5 or 6 (depending on f0, f1, or f2) coordinates should be plotted. 
%    
%   Inputs:
%       obj - A RegCRTBPAtlas object  
%       evalNode - {S, T} are both vectors of space and time evaluation nodes respectively. These nodes are computed in the GLOBAL parameterization domain.  
%       chartRegType - An element of {0, 1, 2} to define which charts should be plotted. Use repeated calls for handling multiple regtypes in one plot. 
%       plotRegType - An element of {0, 1, 2} to define which coordinates to plot in. 
%       plotIdx - Use [1,3] for configuration space plots. For plotting one velocity component use [1, 3, 2] or [1, 3, 4].
%
%   Outputs:
%       none
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 06-Feb-2020; Last revision: 21-Apr-2020

%% parse input 
% parse input and varargin
p = inputParser;
addRequired(p, 'obj');
addRequired(p, 'evalNode');
addRequired(p, 'chartRegType');
addRequired(p, 'plotRegType');
addRequired(p, 'plotIdx');
addParameter(p, 'FaceColor', 'interp')
addParameter(p, 'EdgeColor', 'none')
addParameter(p, 'FaceAlpha', 1)

parse(p, obj, evalNode, chartRegType, plotRegType, plotIdx, varargin{:})
faceColor = p.Results.FaceColor;
edgeColor = p.Results.EdgeColor;
faceAlpha = p.Results.FaceAlpha;

%% set up patches for atlas
globalSpace = evalNode{1};
globalTime = evalNode{2};
atlasDimension = length(evalNode);
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

% filter out charts of the correct regtype
regIdx = arrayfun(@(chart)isequal(chart.RegType, chartRegType), obj.Chart); % logical indices for charts in correct regtype
regChart = obj.Chart(regIdx); % slice all charts in the correct regtype

% define change of coordinate map which maps points in chartReg to points in the plot regtype which
% may be different.
if isequal(chartRegType, plotRegType)
    map2plotreg = @(x)x; % identity map
else
    map2plotreg = @(x)CRTBP2reg(x(:,1:4), obj.Chart(1).Parameter(1), plotRegType - chartRegType);
end

clear iPatchCell iPatchArray iPatch2PlotRegType iPatchEval
% main loop
for iChart = regChart % loop through charts of correct regtype
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
        iPatchArray = cell2mat(iPatchCell); % n-by-k where n = 5 (f0) or 6 (f1 and f2)
        iPatch2PlotRegType = map2plotreg(iPatchArray); % map each evaluation point to plot regtype
        iPatchEval = iPatch2PlotRegType(:, plotIdx); % slice out coordinates to plot
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
if strcmp(faceColor, 'interp')
    patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData, 'FaceColor', 'interp', 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
else
    patch('Faces', face, 'Vertices', vertex, 'FaceColor', faceColor, 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
end

% patch('Faces', face, 'Vertices', vertex, 'FaceColor', patchColor, 'EdgeColor', 'none', 'FaceAlpha', .5);
% patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData, 'FaceColor', 'interp', 'EdgeColor', 'none');
% patch('Faces', face, 'Vertices', vertex, 'FaceColor',patchColor,'EdgeColor',edgeColor,'FaceAlpha',.7);
end % end patch

% Revision History:
%{
21-Apr-2020: Renamed to patch from CRTBPpatch. This is now an overloaded version of the patch method for the RegCRTBPAtlas class
%}
