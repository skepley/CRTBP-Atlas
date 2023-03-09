function t = taylorregtime(obj, varargin)
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
%   email: s.kepley@vu.nl
%   Date: 22-Apr-2020;

%% parse input

if isequal(nargout, 0)
    warning('Check the calling function. This function was called with nargout = 0. This used to be called to update the regTime for the chart but this should no longer be used.')
end

if nargin > 1
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
        if nargin > 2 % compute integrand and antiderivative as a polynomial and then evaluate
            error('This formula may not be trustworthy')
            p = fliplr(conv(x_i, x_i).');
            q = fliplr(conv(y_i, y_i).');
            t = 4*tau*(polyval(polyint(p + q), tau/obj.Tau));
            
        else  % integral on [0,1] is just the weighted sum of coefficients
            t = 4*tau*(sum(conv(x_i, x_i)./(1:2*M-1).') + sum(conv(y_i, y_i)./(1:2*M-1).')); % integrate dt/dTau_i
        end
end

end % end taylorregtime

% Revision History:
%{
13-Aug-2020 - Added a version which returns the F0-time instead of updating the chart properties. This occurs if called
with an output argument. Also added an optional argument to specify a F1/F2 time instead of using the charts terminal time.
This is used to evaluate orbits.
20-May-2021 - Fixed the formula for evaluating t(tau) for some tau < obj.Tau
8-March-2023 - Rewritten entirely for simplicity. It should only be used to retrieve regtimees and never to update Chart properties.
%}
