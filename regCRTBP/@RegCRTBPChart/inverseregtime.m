function regTau = inverseregtime(obj, F0Tau)
%INVERSEREGTIME - Invert the taylorregtime function.
%
%   INVERSEREGTIME() maps the time of flightn for an orbit segment w.r.t F0 time to the corresponding F1/F2 time of flight
%
%   Syntax:
%       regTau = INVERSEREGTIME(RegCRTBPChart, F0Tau) returns the timestep in regularized time RegTau for this chart
%           given a specified timestep in F0 time.
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 11-Aug-2020; Last revision: 20-Aug-2020


%% parse input
switch obj.RegType
    case 0 % chart is already in F_0 coordinates
        regTau = F0Tau;
        
    otherwise % chart is in f_1 or f_2
        x_i = obj.Coordinate(1).Coefficient(:, 1); % Taylor expansion in time for the orbit through x_i(0)
        y_i = obj.Coordinate(3).Coefficient(:, 1); % Taylor expansion in time for the orbit through y_i(0)
        M = obj.Truncation(1); % Truncation size in time
        P = zeros(2*M, 1); % the map to invert has the form: F0Tau = P(regTau) where p is a univariate polynomial of degree 2M
        P(2:end) = 4*(conv(x_i, x_i) + conv(y_i, y_i))./((1:2*M-1)'.*obj.LocalTime.^(0:2*M-2)');
        P(1) = -F0Tau;  % Substract F0Tau from the constant term
        P = flip(P); % switch to desecnding order to use polyval function
        
        % Solve P(u) - F0Tau = 0
        initTau = obj.LocalTime*F0Tau/obj.Tau; % initial guess
        F = @(x)polyval(P, x);  % zero finding map
        try
            regTau = fzero(F, initTau);
        catch
            sprintf('Inverse time fzero fail')
        end
        %         dF = @(x)polyval(polyder(P), x);  % derivative of zero finding map
        %         regTau = findroot(F, dF, initTau);  % Newton's method to find solutions
        if isnan(regTau)
            warning('inverse regtime rootfinder failed to converge')
        end
end


end % end inverseregtime

% Revision History:
%{
20-Aug-2020 - Switch rootfinder to fzero to handle multiple roots when calling this function with F0Tau = 0.
%}
