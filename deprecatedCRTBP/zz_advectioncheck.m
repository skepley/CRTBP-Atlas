function [addBoundary, addEvaluation] = advectioncheck(obj, advectionChart, varargin)
%ADVECTIONCHECK - Advection check algorithm for RegCRTBPChart atlases
%
%   ADVECTIONCHECK() - A more detailed description of the function
%
%   Syntax:
%       output = ADVECTIONCHECK(input1, input2)
%       [output1, output2] = ADVECTIONCHECK(input1, input2, input3)
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
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
%   Date: 19-Mar-2019; Last revision: 19-Mar-2019

goodchart = @(chart)true; % placeholder for fitness checker
if nargin > 2
    nBootstrap = varargin{1};
else
    nBootstrap = 0;
end

if ~goodchart(advectionChart)
    % here we describe how we want to subdivide
    addBoundary = subdividerule(advectionChart);
    addEvaluation = {};
else
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
