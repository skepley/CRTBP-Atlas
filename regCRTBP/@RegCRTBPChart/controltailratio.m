function newBd = controltailratio(obj, tailStartsAt, maxTailRatio)
%CONTROLTAILRATIO - For a row vector of boundaryCharts return an array of boundaryCharts parameterizing the same image with given bound on tailratio.
%
%   Syntax:
%       newBoundary = CONTROLTAILRATIO(oldboundary, N, r)  Returns a vector of boundary charts which parameterize the same 
%           image as oldboundary but each chart in the array has at least 100*(1-r)-percent of its ell^1 weight in the first 
%           N-1 coefficients.
%
%   Inputs:
%       tailStartsAt - N denotes which is the first index defined to be the "tail". Smaller yields more conservative subdivision. 
%       maxTailRatio - Ratio in the interval [0, 1]. Smaller means more conservative subdivision.

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 25-Mar-2023


nSubChart = 2;  % number of subdivisions. This exists just in case we want to make it a parameter some day.
boundaryStack = obj; % initialize the stack of unchecked boundary charts.
newBd = RegCRTBPChart; % initialize the output
newBd = newBd(2:end); % empty the output array

while ~isempty(boundaryStack)
    iBoundaryChart = boundaryStack(1); % pop first chart off boundary stack
    tailRatio = iBoundaryChart.Coordinate.tailratio(tailStartsAt);
    if max(tailRatio(~isnan(tailRatio))) > maxTailRatio  % The tailRatio is not controlled. Split and recurse
        boundaryStack = [split_boundary(iBoundaryChart, nSubChart), boundaryStack(2:end)];
        
    else
        newBd = [newBd, iBoundaryChart];  % The tailRatio is already less than the maxTailRatio. Pass to the output array.
        boundaryStack = boundaryStack(2:end);  % remove the passed chart
    end
end
end % end controltailratio

function newBd = split_boundary(bdChart, nSubChart)
endPts = linspace(-1, 1, 1 + nSubChart)';  % partition [-1, 1] into the proper number of uniformly spaced nodes
parmRange = [endPts(1:end-1), endPts(2:end)];  % arrays of pairwise adjacent nodes
newBd = bdChart.subdivide(parmRange);
end
