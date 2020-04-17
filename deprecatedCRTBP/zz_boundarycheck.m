function [addBoundary, addAdvection] = boundarycheck(obj, boundaryChart, maxTau, varargin)
%BOUNDARYCHECK - Boundary check algorithm for RegCRTBPChart atlases
%
%   BOUNDARYCHECK() - Checks the fitness of a boundary chart for the reguarlized CRTBP. The outcome
%   of this check can be one of the following:
%       1. Advection - This chart is well suited to continue to the advection phase.
%       2. Subdivide - This chart violates some condition which can be solved by subdivision.
%       3. Crash - This chart violates some condition which immediately flags it as a crashed chart. This means
%           the chart is a terminal boundary chart for the atlas.
%
%   Syntax:
%       [addBoundary, addAdvection] = BOUNDARYCHECK(atlas, bdChart, maxTau)
%
%   Inputs:
%       obj - An atlas object containing a collection of d-dimensional charts
%       boundaryChart - A (d-1)-dimensional Chart object
%       maxTau - The maximum integration time for which the atlas continues boundaries by advection
%
%   Outputs:
%       addBoundary - A collection of boundary charts to add to the stack. This happens if the boundarychart was subdivided
%       addAdvection - A single boundary chart which has passed all fitness checks. This will continue to the advection step
%           *** ONE OF THESE OUTPUTS SHOULD ALWAYS BE EMPTY ***
%
%   Subfunctions: rule1, rule2, crashnow
%   Classes required: Atlas, Chart, Scalar
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 19-Mar-2019; Last revision: 31-Mar-2019

% subdivisionRule = {@rule1,@rule2}; % cell array of fitness check functions
subdivisionRule = {@rule2}; % cell array of fitness check functions
nRule = length(subdivisionRule);

%% ---------------------- CRASH CONDITIONS ----------------------
% define conditions under which advection of this boundary should immediately stop
if boundaryChart.Crash
    warning('A crashed boundary chart was passed into boundarycheck. This should not happen')
    addAdvection = {}; % pass nothing back out
    addBoundary = {};
    return
end
crashNow = abs(boundaryChart.TimeSpan) >= abs(maxTau); % true if integration time exceeded

if crashNow
    boundaryChart.Crash = true;
    parentChart = boundaryChart.ParentHandle;
    newTau = maxTau - parentChart.TimeSpan(1); % get timestep so that final time is exactly maxTau
    parentChart.scaletime(abs(newTau)); % rescale parentChart to newTau timestep
    disp('max integration time reached')
    obj.CrashStack = [obj.CrashStack, boundaryChart]; % append crashed boundary to the atlas stack
    addBoundary = {}; % return nothing
    addAdvection = {};
else
    %% ---------------------- SUBDIVIDE CONDITIONS ----------------------
    % loop through the list of fitness checks (called rules) to determine if this boundary chart violates any of them
    iRule = 1;
    while length(boundaryChart) == 1 && iRule <= nRule
        boundaryChart = subdivisionRule{iRule}(obj, boundaryChart); % check this subdivision rule
        iRule = iRule + 1;
    end
    
    if isempty(boundaryChart) % boundaryChart failed a crashcheck in some fitness test.
        addBoundary = {}; % pass out nothing
        addAdvection = {};
    elseif length(boundaryChart) > 1 % pass subdivided boundaries back to the boundary stack
        addBoundary = boundaryChart;
        addAdvection = {};
    else % made it through all subdivision checks. pass to the advection stack
        addBoundary = {};
        addAdvection = boundaryChart;
    end
end
end % end boundarycheck


% function crashnow(bdChart)
% %CRASHNOW - Immediately crashes this chart and returns nothing to the boundary or advection stacks
%
% end % crashnow

function newBd = rule1(obj, bdChart)
%RULE1 - Prune portions of the chart which get too close to the singularity

maxSquareRadius = 0.25;
u = bdChart.Coordinate(1);
v = bdChart.Coordinate(3);
evalGrid = linspace(-1,1,20);
squareRadius = mid(u.eval(evalGrid).^2 + v.eval(evalGrid).^2);
if max(squareRadius) > maxSquareRadius % some portion of the chart lies outside the circle
    if any(squareRadius <= maxSquareRadius) % some portion of the chart is still inside the circle
        uu_vv = u*u + v*v; % get Scalar parameterization for u^2 + v^2
        F = mid(uu_vv.Coefficient);
        F(1) = F(1) - maxSquareRadius; % F(x) = u^2 + v^2 - r^2 has roots at the endpoints of intervals where the chart violates the threshold
        
        % get real roots in the interval (-1,1)
        R = roots(flip(F));
        rootIdx = imag(R) == 0 & abs(R) <= 1;
        interiorNodes = sort(R(rootIdx)'); % possible sign changes
        
        % find intervals which violate the condition prune them from the chart
        if isempty(interiorNodes) % no actual violations. This should never actually happen
            warning('chart is certainly not tangent to boundary circle so there must be an error here')
            newBd = bdChart;
            return
        else
            % get endpoints of intervals for which F has a sign change
            endPts = [-1,interiorNodes,1]; % all candidates
            chkPts = .5*(endPts(1:end-1) + endPts(2:end)); % midpoints of each interval
            getSign = sign(polyval(flip(F),chkPts));
            signChange = getSign(2:end) - getSign(1:end-1); % pos ---> neg if diff = -2; neg ---> pos if diff = 2;
            leftPts = interiorNodes(signChange == -2);
            rightPts = interiorNodes(signChange == 2);
            
            % check the sign of first and last subintervals
            if getSign(1) < 0 % negative on 1st subinterval
                leftPts = [-1,leftPts];
            end
            if getSign(end) < 0 % negative on last subinterval
                rightPts = [rightPts,1];
            end
            
            % specify subdivision boundary nodes and then subdivide
            parmRange = [leftPts',rightPts'];
            if length(leftPts) > 3
                fprintf('This idiot wants to subdivide into %d pieces \n',size(parmRange,1))
            else
                fprintf('This one tried to leave the circle \n')
            end
            newBd = bdChart.subdivide(parmRange);
        end
    else % entire chart lies outside circle
        % flag this boundary as a crashed chart and append to the crash stack
        bdChart.Crash = true;
        obj.CrashStack = [obj.CrashStack, bdChart];
        newBd = {};
    end
    
else % this boundary chart is inside the circle
    newBd = bdChart;
end
end % rule1

function newBd = rule2(obj, bdChart)
%RULE2 - Subdivides into smaller charts if the tailratio grows too large

% subdivision parameters
tailStartsAt = 4; % The first coefficient which is considered part of the tail
maxTailRatio = .01; % Threshold value for allowable tailratio
nSubChart = 4; % Number of charts to subivide into
% THIS SUBDIVISION CHOOSES SUBDIVISION NODES UNIFORMLY

tailRatio = bdChart.Coordinate.tailratio(tailStartsAt);
if max(tailRatio(~isnan(tailRatio))) > maxTailRatio
    endPts = linspace(-1,1,1 + nSubChart)';
    parmRange = [endPts(1:end-1),endPts(2:end)];
    newBd = bdChart.subdivide(parmRange);
    fprintf('subdivided \n')
else
    newBd = bdChart;
end
end % rule2

% Revision History:
%{
31-Mar-2019 - Added subdivision rule for charts which leave the optimal circle.
%}
