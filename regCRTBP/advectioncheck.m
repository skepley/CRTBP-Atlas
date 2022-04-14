function [addBoundary, addEvaluation] = advectioncheck(obj, advectionChart, varargin)
%ADVECTIONCHECK - Advection check algorithm for RegCRTBPChart atlases
%
%   ADVECTIONCHECK() - Checks the fitness of an interior LorenzChart after advection. The outcome
%   of this check can be one of the following:
%       1. Next Generation - This chart is well conditioned and will be added to the atlas.
%           Its boundary continues to the next generation boundary phase.
%       2. Subdivide - This chart violates some condition which can be solved by subdivision.
%       3. Crash - This chart violates some condition which immediately flags it as a crashed chart. This means
%           the chart is added to the atlas but its boundary is automatically a terminal boundary chart for the atlas.
%
%   Syntax:
%       [addBoundary, addAdvection] = ADVECTIONCHECK(atlas, bdChart, maxTau)
%
%   Inputs:
%       obj - An atlas object containing a collection of RegCRTBPChart objects
%       advectionChart - A d-dimensional RegCRTBPChart object
%
%   Outputs:
%       addBoundary - A collection of boundary charts to add to the stack. This happens if the interior chart failed a check. Its boundary is
%           sent back to the boundary stack for subdivision.
%       addAdvection - An interior chart which has passed all fitness checks. This will continue to the next generation.
%           *** ONE OF THESE OUTPUTS SHOULD ALWAYS BE EMPTY ***
%
%   Subfunctions: none
%   Classes required: @RegCRTBPChart, @Atlas
%   Other m-files required: none 
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 19-Mar-2019; Last revision: 19-Apr-2019

%% ---------------------- CRASH CONDITIONS ----------------------

%% ---------------------- SUBDIVIDE CONDITIONS ----------------------
goodchart = @(chart)true; % placeholder for fitness checker
if ~goodchart(advectionChart) % interior chart failed subdivision fitness test
    addBoundary = subdividerule(advectionChart);
    addEvaluation = {};
else % chart passes advection fitness tests and goes to evaluation stage
    addBoundary = {};
    addEvaluation = advectionChart;
    if advectionChart.IsValid
        warning('shrink-wrapping is still not implemented yet')
        advectionChart.validate(true);
    end
    obj.Chart = [obj.Chart, advectionChart]; % append advected chart to the atlas of charts
    obj.Size = obj.Size + 1; % update atlas size
end
end % end advectioncheck

% Revision History:
%{

%}
