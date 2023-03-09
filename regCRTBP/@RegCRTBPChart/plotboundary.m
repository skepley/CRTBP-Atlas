function varargout = plotboundary(obj, numNode, plotRegType, plotIdx, varargin)
%PLOTBOUNDARY - Plot the boundary lines for a CRTBP Chart
%
%   PLOTBOUNDARY() - A more detailed description of the function
%
%   Syntax:
%       output = PLOTBOUNDARY(input1, input2)
%       [output1, output2] = PLOTBOUNDARY(input1, input2, input3)
%
%   Inputs:
%       input1 - Description
%       numNode - [nT, nS] specifies the number of nodes in each dimension.
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
%   email: s.kepley@vu.nl
%   Date: 06-Feb-2020; Last revision: 06-Feb-2020

% parse input and varargin
p = inputParser;
addRequired(p,'obj');
addRequired(p,'numNode');
addRequired(p,'plotRegType');
addRequired(p,'plotIdx');
addParameter(p,'PlotOptions',{})

parse(p, obj, numNode, plotRegType, plotIdx, varargin{:})
plotOptions = p.Results.PlotOptions;

% set up boundary domain and evaluate chart
timeNode = linspace(0, 1, numNode(1))';
spaceNode = linspace(-1, 1, numNode(2))';
bottom = [spaceNode, timeNode(1)*ones(size(spaceNode))];
right = [spaceNode(end)*ones(size(timeNode)), timeNode];
top = [flipud(spaceNode), timeNode(end)*ones(size(spaceNode))];
left = [spaceNode(1)*ones(size(timeNode)), flipud(timeNode)];
bdEvalData = [bottom; right; top; left];
bdEvalImage = cell2mat(obj.eval(bdEvalData));

if ~isequal(obj.RegType, plotRegType)
    % map which maps points in chartReg to points in the plot regtype
    bdEvalImage = CRTBP2reg(bdEvalImage(:,1:4), obj.Parameter(1), plotRegType - obj.RegType);
end

% plot the boundary data
plotDimension = length(plotIdx);
plotData = mat2cell(bdEvalImage(:, plotIdx), size(bdEvalImage, 1), ones(1, plotDimension));


if isequal(plotDimension, 2)
    h = plot(plotData{:}, plotOptions{:});
elseif isequal(plotDimension, 3)
    h = plot3(plotData{:}, plotOptions{:});
else
    error('Only works for 2 or 3 dimensional indices')
end

% return plot handle if requested
if nargout > 0
    varargout{1} = h;
end



end % end plotboundary

% Revision History:
%{

%}
