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
%   email: s.kepley@vu.nl
%   Date: 19-Mar-2019; Last revision: 12-Jun-2019


%% parse input
p = inputParser;
addRequired(p, 'obj')
addRequired(p, 'boundaryChart')
addRequired(p, 'maxTau')
addParameter(p, 'SubDivideParameter', [4, .1, 4])
addParameter(p, 'HotSwap', true) % by default the charts will automatically apply the conjugacy maps to improve numerical conditioning


parse(p, obj, boundaryChart, maxTau, varargin{:})
subDivideParameter = p.Results.SubDivideParameter;
hotSwap = p.Results.HotSwap;

% cell array of fitness check functions
if hotSwap
    subdivisionRule = {@(obj, boundaryChart)rule1(obj, boundaryChart, subDivideParameter),...
        @(obj, boundaryChart)rule2(obj, boundaryChart, subDivideParameter)}; % check hotswap and tailratio
else
    subdivisionRule = {@(obj, boundaryChart)rule2(obj, boundaryChart, subDivideParameter)}; % turn off hotswapping
end
nRule = length(subdivisionRule);

%% ---------------------- CRASH CONDITIONS ----------------------
% define conditions under which advection of this boundary should immediately stop

% bdChart has previously crashed
if boundaryChart.Crash
    warning('A crashed boundary chart was passed into boundarycheck. This should not happen')
    addAdvection = {}; % pass nothing back out
    addBoundary = {};
    return
end

% bdChart has exceeded the integration time limit
crashNow = abs(boundaryChart.TimeSpan) >= abs(maxTau) || boundaryChart.Generation > obj.MaxGeneration; % true if integration time exceeded or maximum generation exceeded
if crashNow
    parentChart = boundaryChart.ParentHandle;
    
    if abs(boundaryChart.TimeSpan) >= abs(maxTau)  % chart crashed due to maximum integration time
        if isequal(parentChart.RegType, 0)  % local time and regularization time are identical so just rescale time
            newTau = maxTau - parentChart.TimeSpan(1); % get timestep so that final time is exactly maxTau
            parentChart.scaletime(newTau); % rescale parentChart to newTau timestep.
            % 5 July 2019 - newTau updated to have the same sign as obj.Tau so time orientation is preserved correctly.
            
        else  % compute the timestep in regularized time before rescaling time
            newTau = maxTau - parentChart.TimeSpan(1); % get timestep in f0 time
            newRegTime = inverseregtime(parentChart, newTau); % compute corresponding timestep in F1/F2 time
            
            % This block is a hacky fix for the following bug. The Chart.scaletime method operates on the Chart's Tau
            % property in order to access the value of the time rescaling which reduces to last coefficient to the correct
            % norm. But for F1/F2 charts this information is moved to the LocalTime property when taylorregtime is called
            % and the Tau property instead keeps track of the timestep in F0 time. To fix this we reset the Tau Property
            % so that scaletime will act properly, then the final call to taylorregtime below will update both correctly.
            parentChart.Tau = parentChart.LocalTime;
            parentChart.scaletime(newRegTime);  % rescale time so that max time is exactly maxTau
            
        end
        taylorregtime(parentChart); % update local and global time 
%         disp([parentChart.RegType, parentChart.Tau, parentChart.LocalTime])
        
    else
        disp('maximum generation reached')
    end
    
    boundaryChart.Crash = true;
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
end % end if
end % end boundarycheck


% function crashnow(bdChart)
% %CRASHNOW - Immediately crashes this chart and returns nothing to the boundary or advection stacks
%
% end % crashnow

function newBd = rule1(obj, bdChart, subDivideParameter)
%RULE1 - Swap coordinates if this chart crosses its ideal boundary with respect to its regularization type

if ismember('Chebyshev', bdChart.Basis) % regularized coordinates for Chebyshev require sqrt implementation
    newBd = bdChart;
    return
end

switch bdChart.RegType
    case 0 % check if boundary should change to F1 or F2 coordinates
        x0 = bdChart.Coordinate(1);
        y0 = bdChart.Coordinate(3);
        nGridNode = 20;
        evalGrid = linspace(-1,1,nGridNode);
        X0 = mid(x0.eval(evalGrid)); % evaluate at the nodes
        Y0 = mid(y0.eval(evalGrid));
        
        mu = bdChart.Parameter(1);
        S1Check = (X0 - mu).^2 + Y0.^2;
        S2Check = (X0 - mu + 1).^2 + Y0.^2;
        
        % check both ideal domains
        idealRadiusSquare = 1/16;
        if sum(S1Check < idealRadiusSquare) > 0.5*nGridNode % too much of the arc lies in the ideal domain for F1
            newBd = bdChart.hotswap(1); % swap from F0 to F1 coordinates
            newBd = swapcondition(bdChart, newBd, subDivideParameter); % check swap conditioning
            
        elseif sum(S2Check < idealRadiusSquare) > 0.5*nGridNode % too much of the arc lies in the ideal domain for F2
            newBd = bdChart.hotswap(2); % swap from F0 to F2 coordinates
            newBd = swapcondition(bdChart, newBd, subDivideParameter); % check swap conditioning
            
        else % stay in F0 coordinates
            newBd = bdChart;
        end
        
    case 1 % check if boundary should change to F0 coordinates
        idealRadiusSquare = 1/4;
        x1 = bdChart.Coordinate(1);
        y1 = bdChart.Coordinate(3);
        nGridNode = 20;
        evalGrid = linspace(-1,1,nGridNode);
        squareRadius = mid(x1.eval(evalGrid).^2 + y1.eval(evalGrid).^2);
        if sum(squareRadius > idealRadiusSquare) > 0.5*nGridNode % more than half of the arc lies outside the ideal domain
            newBd = bdChart.hotswap(0); % swap from F1 to F0 coordinates
            newBd = swapcondition(bdChart, newBd, subDivideParameter); % check swap conditioning
        else % stay in F1 coordinates
            newBd = bdChart;
        end
        
    case 2 % check if boundary should change to F0 coordinates
        idealRadiusSquare = 1/4;
        x2 = bdChart.Coordinate(1);
        y2 = bdChart.Coordinate(3);
        nGridNode = 20;
        evalGrid = linspace(-1,1,nGridNode);
        squareRadius = mid(x2.eval(evalGrid).^2 + y2.eval(evalGrid).^2);
        if sum(squareRadius > idealRadiusSquare) > 0.5*nGridNode % more than half of the arc lies outside the ideal domain
            newBd = bdChart.hotswap(0); % swap from F2 to F0 coordinates
            newBd = swapcondition(bdChart, newBd, subDivideParameter); % check swap conditioning
        else % stay in F2 coordinates
            newBd = bdChart;
        end
end
end % rule1

function newBd = rule2(obj, bdChart, subDivideParameter)
%RULE2 - Subdivides into smaller charts if the tailratio grows too large

if ismember('Chebyshev', bdChart.Basis) % regularized coordinates for Chebyshev require sqrt implementation
    newBd = bdChart;
    return
end

% unpack subdivision parameters
tailStartsAt = subDivideParameter(1); % The first coefficient which is considered part of the tail
maxTailRatio = subDivideParameter(2); % Threshold value for allowable tailratio
nSubChart = subDivideParameter(3); % Number of charts to subivide into
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

function newBd = swapcondition(oldBd, newBd, subDivideParameter)
%SWAPCONDITIONING - Check if the boundary is well conditioned after a hotswap. If not return a pre-swap subdivision.

% check if swapping vector fields caused the arc to become very poorly conditioned. This happens if the sqrt or inverse
% truncates a lot of significant terms. In this case, the original boundary should be subdivided and passed back to the stack.

% unpack subdivision parameters
tailStartsAt = subDivideParameter(1); % The first coefficient which is considered part of the tail
maxTailRatio = subDivideParameter(2); % Threshold value for allowable tailratio
nSubChart = subDivideParameter(3); % Number of charts to subivide into
tailRatio = newBd.Coordinate.tailratio(tailStartsAt);
if max(tailRatio(~isnan(tailRatio))) > maxTailRatio % check tailratio in new coordinates
    endPts = linspace(-1,1,1 + nSubChart)';
    parmRange = [endPts(1:end-1),endPts(2:end)];
    newBd = oldBd.subdivide(parmRange); % subdivide arc in old coordinates
end
end

% %% OLD CODE SNIPPET FOR PRUNING ARCS
% maxSquareRadius = 0.25;
% u = bdChart.Coordinate(1);
% v = bdChart.Coordinate(3);
% nGridNode = 20;
% evalGrid = linspace(-1,1,nGridNode);
% squareRadius = mid(u.eval(evalGrid).^2 + v.eval(evalGrid).^2);
% if max(squareRadius) > maxSquareRadius % some portion of the chart lies outside the circle
%
%     if any(squareRadius <= maxSquareRadius) % some portion of the chart is still inside the circle
%         uu_vv = u*u + v*v; % get Scalar parameterization for u^2 + v^2
%         F = mid(uu_vv.Coefficient);
%         F(1) = F(1) - maxSquareRadius; % F(x) = u^2 + v^2 - r^2 has roots at the endpoints of intervals where the chart violates the threshold
%
%         % get real roots in the interval (-1,1)
%         R = roots(flip(F));
%         rootIdx = imag(R) == 0 & abs(R) <= 1;
%         interiorNodes = sort(R(rootIdx)'); % possible sign changes
%
%         % find intervals which violate the condition prune them from the chart
%         if isempty(interiorNodes) % no actual violations. This should never actually happen
%             warning('chart is certainly not tangent to boundary circle so there must be an error here')
%             newBd = bdChart;
%             return
%         else
%             % get endpoints of intervals for which F has a sign change
%             endPts = [-1,interiorNodes,1]; % all candidates
%             chkPts = .5*(endPts(1:end-1) + endPts(2:end)); % midpoints of each interval
%             getSign = sign(polyval(flip(F),chkPts));
%             signChange = getSign(2:end) - getSign(1:end-1); % pos ---> neg if diff = -2; neg ---> pos if diff = 2;
%             leftPts = interiorNodes(signChange == -2);
%             rightPts = interiorNodes(signChange == 2);
%
%             % check the sign of first and last subintervals
%             if getSign(1) < 0 % negative on 1st subinterval
%                 leftPts = [-1,leftPts];
%             end
%             if getSign(end) < 0 % negative on last subinterval
%                 rightPts = [rightPts,1];
%             end
%
%             % specify subdivision boundary nodes and then subdivide
%             parmRange = [leftPts',rightPts'];
%             if length(leftPts) > 3
%                 fprintf('This idiot wants to subdivide into %d pieces \n',size(parmRange,1))
%             else
%                 fprintf('This one tried to leave the circle \n')
%             end
%             newBd = bdChart.subdivide(parmRange);
%         end
%     else % entire chart lies outside circle
%         % flag this boundary as a crashed chart and append to the crash stack
%         bdChart.Crash = true;
%         obj.CrashStack = [obj.CrashStack, bdChart];
%         newBd = {};
%     end
%
% else % this boundary chart is inside the circle
%     newBd = bdChart;
% end


% Revision History:
%{
31-Mar-2019 - Added subdivision rule for charts which leave the optimal circle.
19-May-2019 - Changed rule 1 to do automatic coordinate switching based on ideal domains for the 3 vector fields instead of just pruning.
%}
