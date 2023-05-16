function t = physical_time(obj, tau)
%PHYSICAL_TIME - A version of the taylorregtime function which returns the physical time but does not change the chart's Tau or Timespan properties
%
% Inputs:
% obj - Description of the object (chart)
% tau - Value of Tau
%
% Outputs:
% t - Physical time
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 28-Mar-2023; 

switch obj.RegType
    case 0 % chart is already in F0 coordinates
        t = tau;
    otherwise % chart is in f_1 or f_2
        x_i = obj.Coordinate(1).Coefficient(:,1); % Taylor expansion in time for the orbit through x_i(0)
        y_i = obj.Coordinate(3).Coefficient(:,1); % Taylor expansion in time for the orbit through y_i(0)
        M = obj.Truncation(1); % Truncation size in time
        t = 4*tau*(sum(conv(x_i, x_i)./(1:2*M-1).') + sum(conv(y_i, y_i)./(1:2*M-1).')); % integrate dt/dTau_i
end
end % end physical_time

