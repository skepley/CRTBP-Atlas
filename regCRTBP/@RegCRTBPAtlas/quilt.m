function varargout = quilt(obj)
%QUILT - Plot the domain of an atlas of charts with each chart's global domain as a different color
%
%   QUILT() - A more detailed description of the function
%
%   Syntax:
%       h = QUILT(atlas) produces a plot of the domain for this atlas. 
%       
%       QUILT(atlas) create and dispaly the plot without returning any handle
%    
%   Inputs:
%       obj - An Atlas object
%
%   Outputs:
%       h - A handle to the patch plot
%
%   Subfunctions: none
%   Classes required: Atlas
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 29-Mar-2019; Last revision: 29-Mar-2019

h = figure;
set(h, 'OuterPosition', [560   165   980   860])
hold on
% loop over atlas and create a single box patch for each chart's domain using random colors
for chart = obj.Chart
    s = chart.SpatialSpan;
    t = chart.TimeSpan;
    v = [s(1), t(1); s(2), t(1); s(2), t(2); s(1), t(2); s(1), t(1)];
    fill(v(:,1),v(:,2), rand(1,3))
end

% return handle if necessary
if nargout > 0
    varargout{1} = h;
end

end % end quilt

% Revision History:
%{

%}
