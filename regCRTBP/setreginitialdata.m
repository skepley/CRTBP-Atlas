function initialData = setreginitialdata(mu, regType, N, arc, isValid)
%SETREGINITIALDATA - Set up initial data for the collision manifold of the regularized CRTBP field
%
%   SETREGINITIALDATA() - A more detailed description of the function
%
%   Syntax:
%       output = SETREGINITIALDATA(input1, input2)
%       [output1, output2] = SETREGINITIALDATA(input1, input2, input3)
%    
%   Inputs:
%       mu - A value in (0,0.5] which specifies the mass of the small primary.
%       regType - 1 or 2 determines which primary is desingularized.
%       N - Spatial truncation
%       arc - [arcMid, halfArcLength] defines an arc segment to parameterize
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
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jun-2019; Last revision: 18-Jun-2019

%  If mu <= 0.5 % set up some initial data on a circle in the zero energy level set of the F2 field otherwise the initial data lies on the F1 field
arcMid = arc(1);
halfArcLength = arc(2);

% radius is either 8*mu or 8*(1-mu) for the small/large radius respectively
if isequal(regType, 2)
    IDRadius = sqrt(8*mu);
elseif isequal(regType, 1)
    IDRadius = sqrt(8*(1-mu));
end
expCoeff = IDRadius*exp(1i*arcMid)*((1i*halfArcLength).^(0:N-1)./factorial(0:N-1));
xInitial = zeros(1,N);
pInitial = real(expCoeff);
yInitial = zeros(1,N);
qInitial = imag(expCoeff);
initialData = [xInitial; pInitial; yInitial; qInitial];


if isValid
    initialData = intval(initialData);
end

end % end setreginitialdata

% Revision History:
%{

%}
