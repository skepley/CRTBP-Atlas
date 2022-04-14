function unstableBd = L4_local_unstable_boundary(localDataFileName, numUnstableSegment, truncation, mu, tauGuess)
%L4_LOCAL_UNSTABLE_BOUNDARY - Compute a parameterization for the boundary of the local unstable manifold at L4
%
%   L4_LOCAL_UNSTABLE_BOUNDARY() - A more detailed description of the function
%
%   Syntax:
%       output = L4_LOCAL_UNSTABLE_BOUNDARY(input1, input2)
%       [output1, output2] = L4_LOCAL_UNSTABLE_BOUNDARY(input1, input2, input3)
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
%   email: shane.kepley@rutgers.edu
%   Date: 04-Mar-2022; Last revision: 04-Mar-2022

% first lift local L4 unstable manifold into phase space and get boundary

load(localDataFileName) % get local manifold data and pick off coordinates
XuLocal = mid(Au);
PuLocal = mid(Bu);
YuLocal = mid(Cu);
QuLocal = mid(Du);
RuLocal = mid(Ru);
SuLocal = mid(Su);

% take equally spaced nodes on the circle
theta = linspace(0, 2*pi, numUnstableSegment + 1);
nodes = exp(1i*theta);
globalSpatialSpan = (theta - pi)/pi;

% lift parameterization
clear unstableBd;
for k = 1:numUnstableSegment
    % get linear segment between nodes
    initNode = nodes(k);
    endNode = nodes(k+1);
    thisSegment(1,:) = [.5*(initNode + endNode),.5*(endNode - initNode)];
    thisSegment(2,:) = conj(thisSegment(1,:));
    
    % lift through Ps
    X0 = bivar_polycomp(XuLocal, thisSegment(1,:), thisSegment(2,:));
    P0 = bivar_polycomp(PuLocal, thisSegment(1,:), thisSegment(2,:));
    Y0 = bivar_polycomp(YuLocal, thisSegment(1,:), thisSegment(2,:));
    Q0 = bivar_polycomp(QuLocal, thisSegment(1,:), thisSegment(2,:));
    R0 = bivar_polycomp(RuLocal, thisSegment(1,:), thisSegment(2,:));
    S0 = bivar_polycomp(SuLocal, thisSegment(1,:), thisSegment(2,:));
    
    % append to local parameterization boundary data
    localUnstableBd = real([X0(1:truncation(2)); P0(1:truncation(2)); Y0(1:truncation(2)); Q0(1:truncation(2)); R0(1:truncation(2)); S0(1:truncation(2))]);
    kSpatialSpan = globalSpatialSpan(k:k+1); % get coordinates in interval [-1,1]
    unstableBd(k) =  RegCRTBPChart(localUnstableBd, 'Taylor', 0, truncation, mu, 0, 'InitialScaling', tauGuess, 'boundary', true);
    unstableBd(k).SpatialSpan = kSpatialSpan;
end

end % end L4_local_unstable_boundary

% Revision History:
%{

%}
