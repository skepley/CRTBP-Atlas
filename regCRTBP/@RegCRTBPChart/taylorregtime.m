function taylorregtime(obj)
%TAYLORREGTIME - A version of regtime which uses the Taylor expansion to compute the time rescaling instead of ode45
%
%   TAYLORREGTIME() - Computes the function t(tau) = int_{t0}^{t0 + tau} 4*w*wBar ds
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
%   Date: 22-Apr-2020; Last revision: 22-Apr-2020

%% parse input
switch obj.RegType
    case 0 % chart is already in F0 coordinates
        return
    otherwise % chart is in F1 or F2
        initialCondition = obj.eval([0,0]); % get an initial point on the chart
        F = @(t,x)rk45regvectorfield(t, x, obj.Parameter, obj.RegType);
        [tau, Fsol] = ode45(F, obj.TimeSpan, [initialCondition{:}]', odeset('AbsTol',1e-13));
        w = Fsol(:,1).^2 + Fsol(:,3).^2;
        F0Tau = 4*trapz(tau,w);
        obj.TimeSpan = [obj.TimeSpan(1), obj.TimeSpan(1) + F0Tau];
        obj.Tau = diff(obj.TimeSpan);
end 
end % end taylorregtime

% Revision History:
%{

%}
