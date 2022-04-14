function varargout = rk45regtime(obj)
%REGTIME - Converts time in F1 or F2 to the equivalent orbit time in F0
%
%   REGTIME() - Computes the function t(tau) = int_{t0}^{t0 + tau} 4*w*wBar ds
%
%   Inputs:
%       obj - An interior chart parameterizing a trajectory for Fi on a time interval [0, tau]
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
%   Date: 21-May-2019; Last revision: 13-Aug-2020

%% parse input
warning('This function has not been updated to allow varargin. It is only for testing taylorregtime is working properly.') 
obj.LocalTime = obj.Tau;  % keep track of local time in order to map back when needed

switch obj.RegType
    case 0 % chart is already in F0 coordinates
        F0Tau = obj.Tau;
        
    otherwise % chart is in F1 or F2
        initialCondition = cell2mat(obj.eval(zeros(1, obj.Dimension(1)))); % get an initial point on the chart
        F = @(t,x)rk45regvectorfield(t, x, obj.Parameter, obj.RegType);
        [tau, Fsol] = ode45(F, obj.TimeSpan, initialCondition.', odeset('AbsTol', 1e-13));
        w = Fsol(:,1).^2 + Fsol(:,3).^2;
        F0Tau = 4*trapz(tau,w);
end

if nargout > 0
    varargout{1} = F0Tau;  % just return the regularized time without updating the chart properties
 
elseif ~isequal(obj.RegType, 0)  % update chart properties with the global time for regularized charts
    obj.TimeSpan = [obj.TimeSpan(1), obj.TimeSpan(1) + F0Tau];
    obj.Tau = diff(obj.TimeSpan);
end
end % end regtime

% Revision History:
%{
13-Aug-2020 - Added a version which returns the F0time instead of updating the chart properties. This occurs if called 
with an output argument. 
%}
