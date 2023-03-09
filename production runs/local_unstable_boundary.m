function [unstableBd, evalBd] = local_unstable_boundary(mu, C, source, target, truncation, tauGuess, varargin)
%LOCAL_UNSTABLE_BOUNDARY - Return a parameterization of a local unstable manifold for L4 or a collision
%
%   LOCAL_UNSTABLE_BOUNDARY() - A more detailed description of the function
%
%   Syntax:
%       output = LOCAL_UNSTABLE_BOUNDARY(input1, input2)
%       [output1, output2] = LOCAL_UNSTABLE_BOUNDARY(input1, input2, input3)
%
%   Inputs:
%       mu - mass ratio parameter
%       C - Energy of the local manifold. This should be constant for the entire manifold
%       source - Either 'P1', 'P2', or a filename where local L4 unstable manifold data can be loaded
%       target - This needs to be specified to determine whether the local manifold should be double covered (in the P1 or P2 case)
%       truncation - [M, N]
%       tauGuess - An initial guess for the first integration step. 0.1 works fairly well here.
%
%   Outputs:
%       unstableBd - A vector of RegCRTBP Boundary charts which parameterize the boundary of the collision or L4 unstable manifold.
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
if any(strcmp(source, {'P1', 'P2'})) && ~any(strcmp(target, {'P1', 'P2'})) % source is a collision and target is L4. This requires double covering the local unstable manifold.
    initialUnstableArc = [pi/2, pi]; % double cover manifold i.e. parameterize the entire local collision manifold (circle) centered at pi/2
    
elseif any(strcmp(source, {'P1', 'P2'}))   % source is a collision but only needs to be single covered
    initialUnstableArc = [pi/2, pi/2]; % single cover manifold i.e. parameterize half the local collision manifold (circle) centered at pi/2
    
else  % source is an L4 unstable manifold
    numUnstableSegment = 8; % Default choice of 8 pieces for L4 local boundary
end


switch source
    case 'P1'
        initialData = setreginitialdata(mu, 1, truncation(2), initialUnstableArc, false); % parameterize the local collision manifold in f_1 coordinates
        unstableBd = RegCRTBPChart(initialData, 'Taylor', 0, truncation, [mu, C], 1, 'InitialScaling', tauGuess, 'boundary', true);  % prepare boundary chart to integrate forward in time
        evalBd = @(s)local_collision_unstable(initialUnstableArc, s);
        
    case 'P2'
        initialData = setreginitialdata(mu, 2, truncation(2), initialUnstableArc, false); % parameterize the local collision manifold in f_2 coordinates
        unstableBd = RegCRTBPChart(initialData, 'Taylor', 0, truncation, [mu, C], 2, 'InitialScaling', tauGuess, 'boundary', true);  % prepare boundary chart to integrate forward in time
        evalBd = @(s)local_collision_unstable(initialUnstableArc, s);
        
    otherwise % Target is L4 which must be loaded from specified file, lifted into phase space, and evaluated along a piecewise polygonal boundary
        load(source) % load local manifold data and pick off coordinates
        XuLocal = mid(Au);
        PuLocal = mid(Bu);
        YuLocal = mid(Cu);
        QuLocal = mid(Du);
        RuLocal = mid(Ru);
        SuLocal = mid(Su);
        
        % take equally spaced nodes on the circle
        theta = linspace(0, 2*pi, numUnstableSegment + 1);  % partition of [0, 2pi]
        nodes = exp(1i*theta);
        globalSpatialSpan = (theta - pi)/pi;  % partition of [-1, 1]
        
        % lift parameterization
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
        evalBd = @(s)local_L4_unstable(numUnstableSegment, s);
        
end
end % end local_unstable_boundary

% Define local unstable manifold boundary evaluation maps and return these with parameters curried.
function localParmCoordinates = local_collision_unstable(unstableArc, globalSpace)
% Local manifold coordinate mapping a global spatial coordinate to its preimage under the
% local parameterization of the unstable collision manifold to P1 or P2.

% unstableArc has the form: [midPoint, arcHalfRadius] for an initial arc segment on P1 collision manifold.
localUnstableParmSpace = (globalSpace + 1)/2;  % distance to globalSpace along local parameterization
localParmCoordinates = localUnstableParmSpace*sum(unstableArc) - (1 - localUnstableParmSpace)*diff(unstableArc);
% angle of the collision with respect to the the local collision manifold. Subtract to reverse the diff order
end

function localParmCoordinates = local_L4_unstable(nSegment, globalSpace)
% Local manifold coordinate mapping a global spatial coordinate to its preimage under the
% local parameterization of the unstable manifold to L4.
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
    localParmCoordinates = [1, 1];1q
    
else
    theta = linspace(0, 2*pi, nSegment + 1);
    nodes = exp(1i*theta);  % equally spaced nodes on the circle
    globalSpatialSpan = (theta - pi)/pi;  % nodes projected to the interval [-1, 1]
    rightNodeIdx = find(globalSpatialSpan > globalSpace, 1);  % return first index which exceeds the global spatial connection coordinate
    leftNodeIdx = rightNodeIdx-1;
    localUnstableParmSpace = (globalSpace - globalSpatialSpan(leftNodeIdx))/(globalSpatialSpan(rightNodeIdx) - globalSpatialSpan(leftNodeIdx));
    % distance along one segment of the polygonal boundary where the collision occurs
    unstableSigma = (1 - localUnstableParmSpace).*nodes(leftNodeIdx) + localUnstableParmSpace.*nodes(rightNodeIdx);
    localParmCoordinates = [unstableSigma, conj(unstableSigma)];
end
end