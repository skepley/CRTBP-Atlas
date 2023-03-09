%FINAL1 - Production run for connections between L4 and P1 at equal masses.
%   Connecting orbits for the following example:
%       Parameters: equal mass
%       Energy: 3
%       Source: L4
%       Target: P1
%       Time: +5, -5 units
%
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 30-Dec-2019; Last revision: 3-Feb-2021

%% ================================================== SET DATA PATHS AND FILENAMES  ==================================================
clear vars
close all
clc

% set paths for accessing and saving data
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

% set file names to save manifold and collision data
runID = 't4_CL4_4021_equalmass_default'; % t = 5, C = 3, (M,N) = (40,20), mu = 1/2, subdivision = (4,.1, 4)
unstableFileID = ['L4_unstable_', runID]; % file or folder name to save unstable manifold atlas data
stableFileID = ['P1_stable_', runID]; % file or folder name to save stable manifold atlas data
collisionFileID = [savePath, 'L4_P1_connections_', runID, '.mat']; % filename to save collision data

%% ================================================== SET RUN PARAMETERS  ==================================================
% Manifold run parameters
mu = 0.5; % equal masses means m1 = m2 = 1/2
C = 3; % L4 energy for equal masses
minTau = -4; % stable manifold integration time
maxTau = 4; % unstable manfiold integration time
numUnstableSegment = 8; % number of initial subdivisions for unstable boundary
initialStableArc = [pi/2, pi]; % [midPoint, arcRadius] for initial on P1 collision manifold. [pi/2, pi/2] parameterizes only a half circle
%   to take advantage of antipodal symmetry for equal masses
regTime = @taylorregtime;  % set method of integration for regularizing time

% Integrator parameters
N = 21; % spatial truncation
M = 40; % temporal truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
tauGuess = 0.1; % initial timestep to scale each Taylor step

% local stable/unstable preimage mappings (curry run parameters)
unstableLocalMap = @(s)localunstablecoordinates(numUnstableSegment, s);
stableLocalMap = @(s)localstablecoordinates(initialStableArc, s);

%%  ================================================== GROW OR LOAD UNSTABLE MANIFOLD FOR L4  ==================================================
if exist([savePath, unstableFileID]) || exist([savePath, unstableFileID, '.mat']) % check if unstable manifold is computed already
    disp('Using existing unstable manifold data')
    atlasL4F = chunkyload(savePath, unstableFileID);
    
else % compute and save unstable manifold data
    disp('Computing new unstable manifold data')
    unstableBd = L4_local_unstable_boundary('newL4_equalmass_local_manifolds', numUnstableSegment, truncation, mu, tauGuess);
    
    % Integrate the unstable manifold forward in time.
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasL4F = RegCRTBPAtlas(unstableBd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasL4F.LeafStack)
        %         fprintf('%0.4f \n', min([atlasL4F.LeafStack.TimeSpan]))
        atlasL4F.growboundary('RegTime', regTime)
    end
    chunkysave(atlasL4F, savePath, unstableFileID) % save manifold data
end


%%  ================================================== GROW OR LOAD STABLE MANIFOLD FOR P1  ==================================================
if exist([savePath, stableFileID]) || exist([savePath, stableFileID, '.mat'])% check if unstable manifold is computed already
    disp('Using existing stable manifold data')
    atlasP1B = chunkyload(savePath, stableFileID);
    
else % compute and save stable manifold data
    disp('Computing new stable manifold data')
    initialData = setreginitialdata(mu, 1, N, initialStableArc, false); % set up some initial data on a F1 circle
    
    % Integrate the collision manifold for the large primary backwards in time.
    stableBd = RegCRTBPChart(initialData, 'Taylor', 0, truncation, [mu, C], 1, 'InitialScaling', -tauGuess, 'boundary', true);
    bdCheck = @(obj, boundaryChart, minTau)boundarycheck(obj, boundaryChart, minTau, 'SubDivideParameter', subDivideParameter);
    atlasP1B = RegCRTBPAtlas(stableBd, -tauGuess, bdCheck, @advectioncheck, 'MaxTau', minTau);
    while ~isempty(atlasP1B.LeafStack)
        atlasP1B.growboundary('RegTime', regTime)
    end
    chunkysave(atlasP1B, savePath, stableFileID) % save manifold data
end

%% ================================================== MINE FOR INTERSECTIONS ==================================================
if isfile(collisionFileID) % check if collisions have already been computed
    disp('Using existing collision data')
    load(collisionFileID)
    
else
    % load stable/unstable manifolds
    atlasB = atlasP1B; % set (stable) backward time atlas variable
    atlasF = atlasL4F; % set (unstable) forward time atlas variable
    sols = mine_intersections(atlasB, atlasF);
    save(collisionFileID, 'sols') % save collision data
end

%% ================================================== LOCAL FUNCTIONS ==================================================

% define coordinate maps for preimages of the local stable/unstable manifolds to map connections back to the BVP parameter
function localParmCoordinates = localunstablecoordinates(nSegment, globalSpace)
% LOCALUNSTABLECOORDINATES Local manifold coordinate mapping a global spatial coordinate to its preimage under the
% local parameterization of the unstable manifold to L4.

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

function localParmCoordinates = localstablecoordinates(stableArc, globalSpace)
% LOCALSTABLECOORDINATES Local manifold coordinate mapping a global spatial coordinate to its preimage under the
% local parameterization of the stable collision manifold to P1.

% stableArc has the form: [midPoint, arcHalfRadius] for an initial arc segment on P1 collision manifold.
localStableParmSpace = (globalSpace + 1)/2;  % distance to globalSpace along local parameterization
localParmCoordinates = localStableParmSpace*sum(stableArc) - (1 - localStableParmSpace)*diff(stableArc);
% angle of the collision with respect to the the local collision manifold. Subtract to reverse the diff order
end