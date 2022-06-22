function [stableBd, evalBd] = local_stable_boundary(mu, C, target, truncation, tauGuess, varargin)
%LOCAL_STABLE_BOUNDARY - Return a parameterization of a local stable manifold for L4 or a collision
%
%   LOCAL_STABLE_BOUNDARY() - A more detailed description of the function
%
%   Syntax:
%       output = LOCAL_STABLE_BOUNDARY(input1, input2)
%       [output1, output2] = LOCAL_STABLE_BOUNDARY(input1, input2, input3)
%
%   Inputs:
%       mu - mass ratio parameter
%       C - Energy of the local manifold. This should be constant for the entire manifold
%       target - Either 'P1', 'P2', or a filename where local L4 stable manifold data can be loaded
%       truncation - [M, N]
%       tauGuess - An initial guess for the first integration step. 0.1 works fairly well here.
%
%   Outputs:
%       stableBd - A vector of RegCRTBP Boundary charts which parameterize the boundary of the collision or L4 stable manifold.
%           these are ready to be integrated backwards in time using the Atlas methods.
%       evalBd - A function which maps from a global space parameter, s in [-1, 1] into the local boundary coordinates. The latter is either
%           a point on a polygonal boundary for L4 or a point on the circle defined by collisions for P1 or P2.
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 03-Jun-2022;

% choose parameterization properties depending on the type of manifold
if any(strcmp(target, {'P1', 'P2'}))  % target is a collision manifold
    initialStableArc = [pi/2, pi]; % double cover manifold i.e. parameterize the entire local collision manifold (circle) centered at pi/2
else  % target is an L4 stable manifold
    numStableSegment = 8; % Default choice of 8 pieces for L4 local boundary
end


switch target
    case 'P1'
        initialData = setreginitialdata(mu, 1, truncation(2), initialStableArc, false); % parameterize the local collision manifold in f_1 coordinates
        stableBd = RegCRTBPChart(initialData, 'Taylor', 0, truncation, [mu, C], 1, 'InitialScaling', tauGuess, 'boundary', true);  % prepare boundary chart to integrate backward in time
        evalBd = @(s)local_collision_stable(initialStableArc, s);
        
    case 'P2'
        initialData = setreginitialdata(mu, 2, truncation(2), initialStableArc, false); % parameterize the local collision manifold in f_1 coordinates
        stableBd = RegCRTBPChart(initialData, 'Taylor', 0, truncation, [mu, C], 2, 'InitialScaling', tauGuess, 'boundary', true);  % prepare boundary chart to integrate backward in time
        evalBd = @(s)local_collision_stable(initialStableArc, s);
        
    otherwise % Target is L4 which must be loaded from specified file, lifted into phase space, and evaluated along a piecewise polygonal boundary
        load(target) % load local manifold data and pick off coordinates
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
            stableBd(k) =  RegCRTBPChart(localStableBd, 'Taylor', 0, truncation, mu, 0, 'InitialScaling', tauGuess, 'boundary', true);
            stableBd(k).SpatialSpan = kSpatialSpan;
        end
        evalBd = @(s)local_L4_stable(numStableSegment, s);
        
end
end % end local_stable_boundary

% Define local stable manifold boundary evaluation maps and return these with parameters curried.
function localParmCoordinates = local_collision_stable(stableArc, globalSpace)
% Local manifold coordinate mapping a global spatial coordinate to its preimage under the
% local parameterization of the stable collision manifold to P1 or P2.

% stableArc has the form: [midPoint, arcHalfRadius] for an initial arc segment on P1 collision manifold.
localStableParmSpace = (globalSpace + 1)/2;  % distance to globalSpace along local parameterization
localParmCoordinates = localStableParmSpace*sum(stableArc) - (1 - localStableParmSpace)*diff(stableArc);
% angle of the collision with respect to the the local collision manifold. Subtract to reverse the diff order
end

function localParmCoordinates = local_L4_stable(nSegment, globalSpace)
% Local manifold coordinate mapping a global spatial coordinate to its preimage under the
% local parameterization of the stable manifold to L4.
%
%Example: For nSegment = 8 this parameterization wraps the segment [-1, 1] around the complex unit circle. In particular:
% P(-1) = 1, P(-0.75) = 0.7071 + 1i*0.7071, P(-0.5) = 1i, P(-0.25) = 0.7071 - 1i*0.7071,
% P(0) = -1, P(0.25) = -0.7071 - 1i*0.7071, P(0.5) = -1i, P(0.75) = 0.7071 - 1i*0.7071.
% The case P(1) has to be handled as an edge case.
%
% visualize with z = arrayfun(P, linspace(-1, 1, 100), 'UniformOutput', false)
%       zz = cell2mat(z)
%       scatter(real(zz), imag(zz))

if isequal(globalSpace, 1)
    localParmCoordinates = [1, 1];
    
else
    theta = linspace(0, 2*pi, nSegment + 1);
    nodes = exp(1i*theta);  % equally spaced nodes on the circle
    globalSpatialSpan = (theta - pi)/pi;  % nodes projected to the interval [-1, 1]
    rightNodeIdx = find(globalSpatialSpan > globalSpace, 1);  % return first index which exceeds the global spatial connection coordinate
    leftNodeIdx = rightNodeIdx-1;
    localStableParmSpace = (globalSpace - globalSpatialSpan(leftNodeIdx))/(globalSpatialSpan(rightNodeIdx) - globalSpatialSpan(leftNodeIdx));
    % distance along one segment of the polygonal boundary where the collision occurs
    stableSigma = (1 - localStableParmSpace).*nodes(leftNodeIdx) + localStableParmSpace.*nodes(rightNodeIdx);
    localParmCoordinates = [stableSigma, conj(stableSigma)];
end
end

