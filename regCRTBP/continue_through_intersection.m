function [sData, uData] = continue_through_intersection(connectionData)
%CONTINUE_THROUGH_INTERSECTION - Given a candidate connection, continue evaluation in time of the stable/unstable charts 
% where the intersection was found. When plotted against the full orbit this should be tangent for a true connection and transverse
% for pseudo-connections
%
%   Syntax:
%       output = CONTINUE_THROUGH_INTERSECTION(input)
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
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
%   Date: 24-Jan-2023; 



stableChart = connectionData{1};
unstableChart = connectionData{2};
localIntersection = connectionData{3};  % x = (s1, t1, s2) and S(s1, t1) = U(s2, t0) in 3 coordinates
localStableSpace = localIntersection(1); % local spatial coordinate for the stable chart at the intersection
localStableRegTau = localIntersection(2); % local time for the stable chart at the intersection
localStableF0Tau = taylorregtime(stableChart, localStableRegTau);
localUnstableSpace = localIntersection(3); % local spatial coordinate for the unstable chart at the intersection
localUnstableF0Tau = 0; % regularize local time with t0 = 0.
globalStableCoords = stableChart.local2global([localStableSpace, localStableF0Tau]);  % global intersection coordinates (Ss, Ts)
globalUnstableCoords = unstableChart.local2global([localUnstableSpace, localUnstableF0Tau]);  % global intersection coordinates (Us, Ut)

localTimeNodes = linspace(0, 1, 100).';
sEvalData = cell2mat(stableChart.eval([globalStableCoords(1)*ones(size(localTimeNodes)), localTimeNodes], 'globalSpace', true));
uEvalData = cell2mat(unstableChart.eval([globalUnstableCoords(1)*ones(size(localTimeNodes)), localTimeNodes], 'globalSpace', true));

sData = sEvalData(:, [1,3]);
uData = uEvalData(:, [1,3]);
end % end continue_through_intersection

