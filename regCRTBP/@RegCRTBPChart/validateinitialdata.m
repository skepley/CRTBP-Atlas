function validateinitialdata(obj)
%VALIDATEINITIALDATA - Appebnd and validate automatic differention terms for initial data.
%
%   Syntax:
%       obj = VALIDATEINITIALDATA(obj) adds the automatica differentiation coordinates to obj and computes rigorous error bounds if the computation is valid.
%
%   Inputs:
%       obj - A RegCRTBPChart with 4 coordinates of initial data specified.
%
%   Outputs:
%       None -
%
%   Subfunctions: none
%   Classes required: RegCRTBPChart
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 19-Mar-2019; Last revision: 10-Apr-2019

% TODO: This needs to be done with validation. The correct way to do it is to solve the ODE: r' = -0.5*u'*r^3 which gives r = 1/sqrt(u).
x = obj.InitialData(1);
y = obj.InitialData(3);
mu = obj.Parameter(1);
mu1 = 1 - mu;

if obj.IsValid
    error('Validation not added yet')
else % just append the automatic differentiation coordinates
    switch obj.RegType
        case 0
            yy = y*y;
            u1 = (x-mu)*(x-mu) + yy;
            R1 = inv(sqrt(u1));
            u2 = (x+mu1)*(x+mu1) + yy;
            R2 = inv(sqrt(u2));
            obj.InitialData(5) = R1;
            obj.InitialData(6) = R2;
        case 1
            xx = x*x;
            yy = y*y;
            u = (xx + yy)*(xx + yy) + 1 + 2*(xx - yy);
            r = inv(sqrt(u));
            obj.InitialData(5) = r;
        case 2
            xx = x*x;
            yy = y*y;
            u = (xx + yy)*(xx + yy) + 1 + 2*(yy - xx);
            r = inv(sqrt(u));
            obj.InitialData(5) = r;
    end
end
end % end validateinitialdata

% Revision History:
%{

%}
