function bestCoordinate = equalmassidealdomain(u, regType)
%EQUALMASSIDEALDOMAIN - Return the ideal domain for a CRTBP position with respect to equal mass parameters.
%
%   EQUALMASSIDEALDOMAIN() - A more detailed description of the function
%
%   Syntax:
%       output = EQUALMASSIDEALDOMAIN(input1, input2)
%       [output1, output2] = EQUALMASSIDEALDOMAIN(input1, input2, input3)
%
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 17-Apr-2020; Last revision: 15-Jun-2020

mu = 1/2;
switch regType
    case 0
        if u(1) < 0
            bestCoordinate = 1;
        else
            bestCoordinate = 2;
        end
        
    case 1
        if 2*u(1).^2 > -mu + sqrt(mu^2 + 4*u(1).^2.*u(2).^2)  % returns true if point is inside the f2 ideal strip
            bestCoordinate = 2;
        else
            bestCoordinate = 1;
        end
        
    case 2
        if 2*u(1).^2 < 1 - mu + sqrt((1-mu)^2 + 4*u(1).^2.*u(2).^2)  % returns true if point is inside the f1 ideal strip
            bestCoordinate = 1;
        else
            bestCoordinate = 2;
        end
end


end % end equalmassidealdomain

% Revision History:
%{

%}
