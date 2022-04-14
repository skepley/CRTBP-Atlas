function growboundary(obj, varargin)
%GROWBOUNDARY - Overloaded method of the growboundary algorithm for the RegCRTBPAtlas class. 
%
%   GROWBOUNDARY() - Performs a single iteration of the following 3 steps: subdivision, advection, and evaluation. This takes k^th generation boundary
%   leaves in the chart tree to the (k+1)^st generation boundary by advecting.
%
%   Syntax:
%       GROWBOUNDARY(obj)
%       GROWBOUNDARY(obj, 'LastCoefficientNorm', 1e-13)
%       GROWBOUNDARY(obj, 'RegTime', func)

%
%   Inputs:
%       obj - A RegCRTBPAtlas instance 
%
%   Outputs: None
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 07-Mar-2019; Last revision: 21-Apr-2020

% parse input
p = inputParser;
addRequired(p, 'obj')
addParameter(p, 'LastCoefficientNorm', eps(1)); % how should the time be rescaled
addParameter(p, 'RegTime', @(chart)chart.regtime()); % global time normalization. By default integration time for the chart is already global. 


parse(p, obj, varargin{:})
lastCoefficientNorm = p.Results.LastCoefficientNorm;
regtime = p.Results.RegTime;

boundarycheck = obj.BoundaryCheck; % inherit boundarychecker
advectioncheck = obj.AdvectionCheck; % inherit advection checker
boundaryStack = obj.LeafStack; % initialize boundary chart stack
advectionStack = obj.ChartClass; % initialize advection stack
advectionStack = advectionStack(2:end); % empty the advection stack
newLeaves = obj.ChartClass; % initialize next generation boundary stack
newLeaves = newLeaves(2:end); % empty the next generation stack 

% main loop
while ~isempty(boundaryStack) % ----------- BOUNDARY PHASE -----------
    iBoundaryChart = boundaryStack(1); % pop first chart off boundary stack
    boundaryStack = boundaryStack(2:end);
    [iBoundary, iAdvection] = boundarycheck(obj, iBoundaryChart, obj.MaxTau); % get new  boundary charts and append to stacks
    boundaryStack = [iBoundary, boundaryStack]; % boundarycheck failures go to front of boundary stack
    advectionStack = [advectionStack, iAdvection]; % boundarycheck pass goes to back of advection stack
    
    while ~isempty(advectionStack) % ----------- ADVECTION PHASE ----------- 
        jAdvectionChart = advectionStack(1); % pop first chart off advection stack
        advectionStack = advectionStack(2:end);
        tau = obj.TimeStepper(jAdvectionChart); % get tau from timestepper
        jAdvectionChart.advect(tau); % convert boundary chart to interior chart by advection
        jAdvectionChart.rescaletime(lastCoefficientNorm); % rescale timestep to force last coefficient below a given norm
        regtime(jAdvectionChart); % normalize time to the global F0 time scale 
        [jBoundary, jEvaluationChart] = advectioncheck(obj, jAdvectionChart); % check if advection chart is good. One output is always empty.
        
        if isempty(jEvaluationChart) % advectioncheck failed
            boundaryStack = [jBoundary, boundaryStack]; % boundary has been subdivided. Return sub-charts to the boundary stack
        else % ----------- EVALUATION PHASE -----------
            newLeaves = [newLeaves, fixtime(jEvaluationChart, 1)]; % evaluate and append to next generation boundary leaves
        end
    end
end
obj.LastGeneration = max([obj.Chart.Generation]);
obj.LeafStack = newLeaves; % replace previous generation boundary stack with next generation boundary stack
fprintf('Atlas Size: %d \n', obj.Size)
fprintf('Generation %d \n', obj.LastGeneration)
fprintf('Next Generation Size: %d \n', numel(obj.LeafStack))
end % end growboundary

% Revision History:
%{

%}
