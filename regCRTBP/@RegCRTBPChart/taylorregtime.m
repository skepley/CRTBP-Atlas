function varargout = taylorregtime(obj, varargin)
%TAYLORREGTIME - A version of regtime which uses the Taylor expansion to compute the time rescaling instead of ode45
%
%   TAYLORREGTIME() - Computes the function t(tau) = int_{t0}^{t0 + tau} 4*w*wBar ds
%
%   Inputs:
%       obj - An interior chart parameterizing a trajectory for F1/F2 on a time interval [0, tau]
%
%   Outputs:
%       F0Time - The equivalent F0 time for this chart
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 22-Apr-2020; Last revision: 13-Aug-2020

%% parse input
if nargin > 1
    if isequal(nargout, 0)
        error('Calling this with a variable argument input should only be done with an output argument.')
    end
    tau = varargin{1};  % F1/F2 time specified as variable argument
else
    tau = obj.Tau;  % compute the F0 time for this chart's terminal time
end

switch obj.RegType
    case 0 % chart is already in F0 coordinates
        t = tau;
    otherwise % chart is in f_1 or f_2
        
        x_i = obj.Coordinate(1).Coefficient(:,1); % Taylor expansion in time for the orbit through x_i(0)
        y_i = obj.Coordinate(3).Coefficient(:,1); % Taylor expansion in time for the orbit through y_i(0)
        M = obj.Truncation(1); % Truncation size in time
        if nargin > 1 % compute integrand and antiderivative as a polynomial and then evaluate
            p = fliplr(conv(x_i, x_i).');
            q = fliplr(conv(y_i, y_i).');
            t = 4*tau*(polyval(polyint(p + q), tau/obj.Tau));
        else  % integral on [0,1] is just the weighted sum of coefficients
            
            t = 4*tau*(sum(conv(x_i, x_i)./(1:2*M-1).') + sum(conv(y_i, y_i)./(1:2*M-1).')); % integrate dt/dTau_i
        end
end

if nargout > 0
    varargout{1} = t;  % just return the regularized time without updating the chart properties
    
else  % update the chart properties with the correct global and local times. This should not be used when tau is
    % specified as a variable argument since it will report the time argument instead of the true terminal chart time.
    obj.LocalTime = tau;  % set the local time to keep track of F1/F2 time interval for this chart.
    if ~isequal(obj.RegType, 0)  % update chart properties with the global time for regularized charts
        obj.Tau = t; % update step size
        obj.TimeSpan = [obj.TimeSpan(1), obj.TimeSpan(1) + t]; % update TimeSpan with correct F0 time
    end
end
end % end taylorregtime

% Revision History:
%{
13-Aug-2020 - Added a version which returns the F0-time instead of updating the chart properties. This occurs if called
with an output argument. Also added an optional argument to specify a F1/F2 time instead of using the charts terminal time.
This is used to evaluate orbits.
20-May-2021 - Fixed the formula for evaluating t(tau) for some tau < obj.Tau
%}
