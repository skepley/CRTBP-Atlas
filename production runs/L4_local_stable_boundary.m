function stableBd = L4_local_stable_boundary(localDataFileName, numStableSegment, truncation, mu, tauGuess)
%L4_LOCAL_STABLE_BOUNDARY - Compute a parameterization for the boundary of the local stable manifold at L4
%
%   L4_LOCAL_STABLE_BOUNDARY() - A more detailed description of the function
%
%   Syntax:
%       output = L4_LOCAL_STABLE_BOUNDARY(input1, input2)
%       [output1, output2] = L4_LOCAL_STABLE_BOUNDARY(input1, input2, input3)
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
%   Date: 06-Mar-2022; Last revision: 06-Mar-2022

% first lift local L4 stable manifold into phase space and get boundary
error('This function was merged into local_stable_boundary.m')
load(localDataFileName) % get local manifold data and pick off coordinates
XsLocal = mid(As);
PsLocal = mid(Bs);
YsLocal = mid(Cs);
QsLocal = mid(Ds);
RsLocal = mid(Rs);
SsLocal = mid(Ss);

% take equally spaced nodes on the circle
theta = linspace(0, 2*pi, numStableSegment + 1);
nodes = exp(1i*theta);
globalSpatialSpan = (theta - pi)/pi;

% lift parameterization
clear stableBd;
for k = 1:numStableSegment
    % get linear segment between nodes
    initNode = nodes(k);
    endNode = nodes(k+1);
    thisSegment(1,:) = [.5*(initNode + endNode),.5*(endNode - initNode)];
    thisSegment(2,:) = conj(thisSegment(1,:));
    
    % lift through Ps
    X0 = bivar_polycomp(XsLocal, thisSegment(1,:), thisSegment(2,:));
    P0 = bivar_polycomp(PsLocal, thisSegment(1,:), thisSegment(2,:));
    Y0 = bivar_polycomp(YsLocal, thisSegment(1,:), thisSegment(2,:));
    Q0 = bivar_polycomp(QsLocal, thisSegment(1,:), thisSegment(2,:));
    R0 = bivar_polycomp(RsLocal, thisSegment(1,:), thisSegment(2,:));
    S0 = bivar_polycomp(SsLocal, thisSegment(1,:), thisSegment(2,:));
    
    % append to local parameterization boundary data
    localStableBd = real([X0(1:truncation(2)); P0(1:truncation(2)); Y0(1:truncation(2)); Q0(1:truncation(2)); R0(1:truncation(2)); S0(1:truncation(2))]);
    kSpatialSpan = globalSpatialSpan(k:k+1); % get coordinates in interval [-1,1]
    stableBd(k) =  RegCRTBPChart(localStableBd, 'Taylor', 0, truncation, mu, 0, 'InitialScaling', -tauGuess, 'boundary', true);
    stableBd(k).SpatialSpan = kSpatialSpan;
end

end % end L4_local_stable_boundary

% Revision History:
%{

%}
