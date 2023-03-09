function sampleData = sample(obj, numNodes)
%SAMPLE - Sample a connecting orbit from the mining algorithm into an approximation for the BVP solver.
%
%   Syntax:
%       output = SAMPLE(input)
%
%   Inputs:
%       obj - RegCRTBPConnection object
%       nNode - Number of uniformly spaced time points to evaluate the orbit. These evaluations are the endpoints of the shooting segments
%
%   Outputs:
%       connection - A struct with the data necessary to form an inital BVP solution guess.
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 23-Feb-2023;


warning('OFF', 'Chart:eval');  % turn off Chart.eval warning messages
warning('OFF', 'Chart:local2global'); 

globalTime = linspace(0, obj.ConnectionTime, numNodes); % a uniform grid of global F0 time points along the entire connecting orbit
globalTime = sort([globalTime, obj.GlobalIntersection(4)]); % Add the intersection time so that it gets evaluated in BOTH the stable and unstable orbit segments.
% This is necessary to fix the problem which occurred when the intersecting charts had different RegTypes.


sampleData = struct; % initialize the connection struct
stableData = obj.sample_connection_from_stable(globalTime);
unstableData = obj.sample_connection_from_unstable(globalTime);
sampleData = stitch_together_connection(sampleData, unstableData, stableData);
end

function sampleData = stitch_together_connection(sampleData, uData, sData)
% Stitch together the stable and unstable segments of a connection and throw out the duplcate intersection if it exists

if isequal(uData.RegVector(end), sData.RegVector(1))  % No hotswap occuring between the intersecting charts. So we merge
    % the last strand of the unstable segment with the first strand of the stable segment.
    
    % trim the redundant data off the stable connection data
    uData.Orbit{end} = cat(2, uData.Orbit{end}, sData.Orbit{1}(:, 2:end));  % throw away duplicate intersection point
    uData.Tau{end} = cat(2, uData.Tau{end}, sData.Tau{1}(2:end)); % throw away redundant zero-time timestep since no hotswap occurred
    sampleData.Orbit = [uData.Orbit, sData.Orbit{2:end}];
    sampleData.Tau = [uData.Tau, sData.Tau{2:end}];
    sampleData.RegVector = [uData.RegVector, sData.RegVector(2:end)];  % throw away redundant regType since the two strands will be merged into one strand
    
else % One of the atlases hotswapped in the generation where this intersection was found
    % We just concatenate all data from the stable/unstable connection evaluations. Wwe have kept BOTH intersection
    % points since they are evaluated on either side of the swap and the timestep between them is zero which is exactly what we want between strands.
    sampleData.Orbit = [uData.Orbit, sData.Orbit];
    sampleData.Tau = [uData.Tau, sData.Tau];
    sampleData.RegVector = [uData.RegVector, sData.RegVector];
end
end % end sample_connection

