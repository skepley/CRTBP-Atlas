function energy = CRTBPenergy(v, parameter, regType)
%ENERGY - One line description of what the function or script performs (H1 line)
%
%   ENERGY() - A more detailed description of the function
%
%   Syntax:
%       output = ENERGY(input1, input2)
%       [output1, output2] = ENERGY(input1, input2, input3)
%
%   Inputs:
%       v - An array whose columns are CRTBP state vectors in R^4
%       parameter - mu or [mu,C] The (small) mass parameter for the CRTBP and the regularization energy if v is in F_1 or F_2 coordinates
%       regType - 0,1 or 2 to specify which potential is to be used for the energy
%
%   Outputs:
%       energy - The value of the energy at this point in phase space. 
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 04-Apr-2019; Last revision: 14-Apr-2019

% unpack coordinates
if ~isequal(size(v,2), 4) % columns of v should be column vectors of points in R^4
    v = v';
end
x = v(:,1);
p = v(:,2);
y = v(:,3);
q = v(:,4);
switch regType
    case 0 % potential for F_0 (standard) coordinates
        mu = parameter(1);
        mu1 = 1-mu;
        r1 = sqrt((x-mu).^2 + y.^2);
        r2 = sqrt((x+mu1).^2 + y.^2);
        U = mu*(0.5*r2.^2 + 1./r2) + mu1*(0.5*r1.^2 + 1./r1);
    case 1 % potential for F_1 (regularized) coordinates
        mu = parameter(1);
        C = parameter(2);
        k1 = x.^2 + y.^2;
        k2 = x.^2 - y.^2;
        r = 1./sqrt(k1.^2 + 1 + 2*k2);
        U = 2*k1.^3 + 4*mu*k1.*k2 + 2*(mu-C).*k1 + 4*(1-mu) + 4*mu*k1.*r;
    case 2 % potential for F_2 (regularized) coordinates
        mu = parameter(1);
        C = parameter(2);
        k1 = x.^2 + y.^2;
        k2 = y.^2 - x.^2;
        r = 1./sqrt(k1.^2 + 1 + 2*k2);
        U = 2*k1.^3 + 4*(mu-1)*k1.*k2 + 2*(1-mu-C).*k1 + 4*mu + 4*(1-mu)*k1.*r;
end
energy = 2*U - p.^2 - q.^2;
end % end energy

% Revision History:
%{

%}
