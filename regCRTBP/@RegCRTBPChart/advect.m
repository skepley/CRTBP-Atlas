function advect(obj, tau)
%ADVECT - Advect a boundary RegCRTBPchart
%
%   Inputs:
%       obj - A boundary type CRTBPchart
%       tau - A timestep to integrate on
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 21-May-2019; Last revision: 21-May-2019

%ADVECT Advect the (d-1) dimensional boundary of a chart to generate the d dimensional interior chart

% Increase Chart dimension and Scalar coordinate dimensions by 1 to prepare for new coefficients to be appended
obj.Dimension(1) = obj.Dimension(1) + 1; % increase domain dimension for Chart

% update Basis with Taylor in time direction
if obj.Dimension > 1
    obj.Basis = {'Taylor', obj.Basis{:}}; 
end

for iCoordinate = 1:obj.Dimension(2) % update dimensions of Scalar coordinates
    obj.Coordinate(iCoordinate).Dimension = obj.Dimension(1);
    obj.Coordinate(iCoordinate).Truncation = [1, obj.Truncation(2:end)];
    obj.Coordinate(iCoordinate).Basis = obj.Basis; 
end
obj.generatecoefficient(tau) % advect this boundary Chart
obj.Tau = tau; % set step size
obj.TimeSpan = [obj.TimeSpan, obj.TimeSpan + tau]; % set time span
end % end advect

% Revision History:
%{
21 May 2019 - Moved out of class definition file and added conversion for F0 time.
20 Feb 2020 - Fixed a bug causing 0 dimensional initial data to fail to compute or compute with 2 terms in obj.Basis 
%}
