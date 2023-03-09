function orbit = condense_bvp_data(bvpData)
%CONDENSE_BVP_DATA - One line description of what the function or script performs (H1 line)
%
%   CONDENSE_BVP_DATA() - A more detailed description of the function
%
%   Syntax:
%       output = CONDENSE_BVP_DATA(input1, input2)
%       [output1, output2] = CONDENSE_BVP_DATA(input1, input2, input3)
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
%   Date: 03-Mar-2022; Last revision: 03-Mar-2022

nSegment = numel(bvpData.Orbit);
orbit = bvpData;
orbit.Tau = zeros(nSegment, 1);
for iSegment = 1:nSegment
    iOrbit = bvpData.Orbit{iSegment};
    orbit.Orbit{iSegment} = iOrbit(:, [1, end]);
    orbit.Tau(iSegment) = sum(bvpData.Tau{iSegment});
end % end condense_bvp_data

% Revision History:
%{

%}
