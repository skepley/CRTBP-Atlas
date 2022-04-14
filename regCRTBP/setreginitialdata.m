function initialData = setreginitialdata(mu, regType, N, arc, isValid)
%SETREGINITIALDATA - Set up initial data for the collision manifold of the regularized CRTBP field
%
%   Syntax:
%       initialData = SETREGINITIALDATA(mu, regType, N, arc, isValid)
%    
%   Inputs:
%       mu - A value in (0,0.5] which specifies the mass of the small primary.
%       regType - 1 or 2 determines which primary is desingularized.
%       N - Spatial truncation
%       arc - [arcMid, halfArcLength] defines an arc segment to parameterize
%       isValid - Specify whether to return validated error bounds for the parameterization (not implemented yet)
%
%   Outputs:
%       initialData - Truncated Taylor expansions for parameterization of the arc segment
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jun-2019; Last revision: 18-Jun-2019

arcMid = arc(1);
halfArcLength = arc(2);

% radius squared is either 8*mu or 8*(1-mu) for the small/large radius respectively
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
