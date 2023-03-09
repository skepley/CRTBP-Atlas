function plot_connection_tangents(obj)
%PLOT_CONNECTION_TANGENTS - Plot the rest of the orbit as tracked through the stable/unstable charts where the intersection was detected
%    
%   Inputs:
%       obj - A connection
%
%   Subfunctions: none
%   Classes required: RegCRTBPConnection
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 01-Mar-2023; 


localTimeNodes = linspace(0, 1, 100).';
sEvalData = cell2mat(obj.StableChart.eval([obj.GlobalIntersection(1)*ones(size(localTimeNodes)), localTimeNodes], 'globalSpace', true));
uEvalData = cell2mat(obj.UnstableChart.eval([obj.GlobalIntersection(3)*ones(size(localTimeNodes)), localTimeNodes], 'globalSpace', true));
sData = sEvalData(:, [1,3]);
uData = uEvalData(:, [1,3]);
plot(sData(:, 1), sData(:, 2), 'r')
plot(uData(:, 1), uData(:, 2), 'k')
end % end plot_connection_tangents

