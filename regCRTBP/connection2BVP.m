function cData = connection2BVP(connectionData, nNode, stableLocalMap, unstableLocalMap)
%CONNECTION2BVP - Process a connecting orbit from the mining algorithm into an approximation for the BVP solver.
%
%   Syntax:
%       output = CONNECTION2BVP(input1, input2)
%       [output1, output2] = CONNECTION2BVP(input1, input2, input3)
%
%   Inputs:
%       connectionData - Cell array of the form: {stableChart, unstableChart, connectionCoordinates}
%           Both Charts are RegCRTBPChart objects which parameterized invariant manifolds and the connection 
%           Coordinates are a vector in R^3 of the form: (s1, t1, s2) such that stable(s1, t1) = unstable(s2, 0)
%       nNode - Number of uniformly spaced time points to evaluate the orbit. These evaluations are the endpoints of the shooting segments
%       stableLocalMap - A function handle which returns the preimage of the local stable manifold parameterization at the endpoint where 
%           it meets the connecting orbit.
%       unstableLocalMap - A function handle which returns the preimage of the local unstable manifold parameterization at the endpoint where 
%           it meets the connecting orbit.
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
%   Date: 20-Aug-2020;

%% ================================================== PARSE CONNECTION PROPERTIES ==================================================
warning('OFF', 'Chart:eval');  % turn off Chart.eval warning messages
warning('OFF', 'Chart:local2global');  % turn off Chart.eval warning messages
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
globalConnectionTime = globalUnstableCoords(2) - globalStableCoords(2);  % time of flight in F0 time
globalTime = linspace(0, globalConnectionTime, nNode); % a uniform grid of global F0 time points along the entire connecting orbit
globalTime = sort([globalTime, globalUnstableCoords(2)]); % Add the intersection time so that it gets evaluated in BOTH the stable and unstable orbit segments.
% This is necessary to fix the problem which occurred when the intersecting charts had different RegTypes.

cData = struct; % initialize the connection struct
cData.ConnectionTime = globalConnectionTime;  % time of flight from one local manifold to the other in F0 time units
cData.UnstableTime = globalUnstableCoords(2);  % time of flight captured in the unstable manifold
% map the local manifold data
cData.LocalUnstable = unstableLocalMap(globalUnstableCoords(1));
cData.LocalStable = stableLocalMap(globalStableCoords(1));

sData = sample_connection_from_stable(stableChart, globalUnstableCoords, globalStableCoords, globalTime);
uData = sample_connection_from_unstable(unstableChart, globalUnstableCoords, globalTime);
cData = stitch_together_connection(cData, uData, sData);

end


function cData = stitch_together_connection(cData, uData, sData)
% Stitch together the stable and unstable segments of a connection and throw out the duplcate intersection if it exists

if isequal(uData.RegVector(end), sData.RegVector(1))  % No hotswap occuring between the intersecting charts. So we merge
    % the last strand of the unstable segment with the first strand of the stable segment.
    
    % trim the redundant data off the stable connection data
    uData.Orbit{end} = cat(2, uData.Orbit{end}, sData.Orbit{1}(:, 2:end));  % throw away duplicate intersection point
    uData.Tau{end} = cat(2, uData.Tau{end}, sData.Tau{1}(2:end)); % throw away redundant zero-time timestep since no hotswap occurred
    cData.Orbit = [uData.Orbit, sData.Orbit{2:end}];
    cData.Tau = [uData.Tau, sData.Tau{2:end}];
    cData.RegVector = [uData.RegVector, sData.RegVector(2:end)];  % throw away redundant regType since the two strands will be merged into one strand
    
else % One of the atlases hotswapped in the generation where this intersection was found
    % We just concatenate all data from the stable/unstable connection evaluations. Wwe have kept BOTH intersection
    % points since they are evaluated on either side of the swap and the timestep between them is zero which is exactly what we want between strands.
    cData.Orbit = [uData.Orbit, sData.Orbit];
    cData.Tau = [uData.Tau, sData.Tau];
    cData.RegVector = [uData.RegVector, sData.RegVector];
end
end
