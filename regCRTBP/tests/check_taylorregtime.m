%CHECK_TAYLORREGTIME - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       CHECK_TAYLORREGTIME description
%
%   Output:
%       CHECK_TAYLORREGTIME output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 22-Apr-2020; Last revision: 22-Apr-2020

clear all
close all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
%============================ MAKE SOME CHOICES ============================
% Integrator parameters
M = 25; % temporal truncation
N = 10; % spatial truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
basis = 'Taylor';
mu = .25; % small mass primary
isValid = false;
hotSwap = false; % specify whether charts should swap based on the ideal domain map or not
bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter,...
    'HotSwap', hotSwap);
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);
%% ============================ 0-DIMENSIONAL DATA, FORWARD TIME ============================
% Integrator time stepping
initialTime = 0;
tauGuess = 0.1; % initial timestep to scale each Taylor step
maxTau = .3;
regType = 1;
% ============================ Compute an Atlas without regularizing time ============================
C = 3;
parameter = [mu,C];

% set up some initial data at a point
p1 = [0.7104, 0.2973, 0.2554, 0.2068];
initialData = p1.';
% Integrate multiple timesteps
ob = RegCRTBPChart(initialData, basis, initialTime, M, parameter, regType, 'InitialScaling', tauGuess);
ob.rescaletime(eps(1));
disp(ob.Tau)

% apply ode45 regtime algorithm
tChart = ob.deepcopy();
tChart.regtime()

% apply Taylor regtime
tauChart = ob.deepcopy();
tauChart.taylorregtime()

disp(tauChart.Tau - tChart.Tau)
return
%% ============================ 1-DIMENSIONAL DATA, FORWARD TIME ============================
% Integrator time stepping
initialTime = 0;
tauGuess = 0.1; % initial timestep to scale each Taylor step
maxTau = .3;
regType = 1;
% ============================ Compute an Atlas without regularizing time ============================
C = 3;
parameter = [mu,C];

% set up some initial data on a line
p1 = [0.7104, 0.2973, 0.2554, 0.2068];
p2 = p1 + 1e-4*ones(size(p1));
initialData = cat(2, 0.5*[p1+p2; p2-p1]', zeros(4,N-2));
% Integrate multiple timesteps
bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
A1 = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
while ~isempty(A1.LeafStack)
    A1.growboundary('RegTime', @(chart)chart.TimeSpan) % do not regularize time
end

% apply Taylor regtime
tauChart = A1.Chart(1).deepcopy();
tauChart.taylorregtime()

% apply ode45 regtime algorithm
tChart = A1.Chart(1).deepcopy();
tChart.regtime()
disp(tauChart.Tau - tChart.Tau)
%% ============================ 1-DIMENSIONAL DATA, BACKWARD TIME ============================

% Integrator time stepping
tauGuess = -0.1; % initial timestep to scale each Taylor step
maxTau = -.3;
regType = 1;

% ============================ Compute an Atlas without regularizing time ============================
C = 3;
parameter = [mu,C];

% set up some initial data on a line
p1 = [0.7104, 0.2973, 0.2554, 0.2068];
p2 = p1 + 1e-4*ones(size(p1));
initialData = cat(2, 0.5*[p1+p2; p2-p1]', zeros(4,N-2));
% Integrate multiple timesteps
bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
A1 = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
while ~isempty(A1.LeafStack)
    A1.growboundary('RegTime', @(chart)chart.TimeSpan) % do not regularize time
end

% apply Taylor regtime
tauChart = A1.Chart(1).deepcopy();
tauChart.taylorregtime()

% apply ode45 regtime algorithm
tChart = A1.Chart(1).deepcopy();
tChart.regtime()
disp(tauChart.Tau - tChart.Tau)


