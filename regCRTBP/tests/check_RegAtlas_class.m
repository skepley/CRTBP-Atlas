%CHECK_REGATLAS_CLASS - A script for developing and testing the RegCRTBPAtlas subclass.
%
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 22-Apr-2020; Last revision: 22-Apr-2020

%% ========================== CHECK SUBCLASS INHERITANCE PROPERTIES AND OVERLOADED METHODS ==========================
clear all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))

% ========================== CHOOSE SOME PARAMETERS ==========================
% Integrator parameters
N = 10; % spatial truncation
M = 25; % temporal truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
basis = 'Taylor';
mu = .25; % small mass primary
isValid = false;
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);

% Integrator time stepping
initialTime = 0;
tDirection = 1;
tauGuess = 0.1; % initial timestep to scale each Taylor step
maxTau = .3;
tf = linspace(0, maxTau, 250); % plotting time points


% ========================== Build an F0 RegCRTBPAtlas ==========================
regType = 0;
parameter = mu;
% set up some initial data on a line
p1 = [0.7104, 0.2973, 0.2554, 0.2068];
p2 = p1 + 1e-4*ones(size(p1));
initialData = cat(2, 0.5*[p1+p2; p2-p1]', zeros(4,N-2));

% set boundary arc
bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);


% Integrate multiple timesteps
bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
A = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
while ~isempty(A.LeafStack)
    A.growboundary()
end

A.patch(3)

%% ================================================== SECTION 2 ==================================================




