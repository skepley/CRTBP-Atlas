function decayEstimate = scalardecaymap(x, M, parameter, regType)
%SCALARDECAYMAP - Evaluates the Taylor coefficient decay map for a 0-dimensional IVP
%
%   Syntax:
%       decayEstimate = SCALARDECAYMAP(x, M, parameter, regType) returns an estimate of the rescaling which makes the
%           last Taylor coefficient have magnitude on order of eps(1).
%
%   Inputs:
%       x - Initial Data
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
%   email: shane.kepley@rutgers.edu
%   Date: 20-Feb-2020; Last revision: 9-Mar-2020

if isequal(regType, 0)
    orbitSegment = RegCRTBPChart(x, 'Taylor', 0, M, parameter, regType);
    orbitSegment.rescaletime(eps(1));
else
    x1Full = CRTBP2reg(reshape(x, [], 4), parameter(1), regType);
    x1 = reshape(x1Full(1:4), [], 1);
    orbitSegment = RegCRTBPChart(x1, 'Taylor', 0, M, parameter, regType);
    orbitSegment.rescaletime(eps(1));
    taylorregtime(orbitSegment);
end
decayEstimate = orbitSegment.Tau;
end % end scalardecaymap

% Revision History:
%{

%}
