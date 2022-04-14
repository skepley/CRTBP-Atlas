%FINAL4 - Production run for ejection-collision connections between P2 and P1 at equal masses for L4 energy.
%   Connecting orbits for the following example:
%       Parameters: equal mass
%       Energy: L4
%       Forward target: P1
%       Backward target: P2
%       Time: +/- 5 units
%
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 02-Jan-2020; Last revision: 02-Jan-2020

%% ================================================== SET DATA PATHS AND FILENAMES  ==================================================
clear all
close all
clc

% set paths for accessing and saving data
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath('/Users/shane/Dropbox/Regularisation3bp/CRTBP Atlas'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

% set file names to save manifold and collision data
runID = 't5_CL4_4020_equalmass_default.mat'; % t = 5, C = 3, (M,N) = (40,20), mu = 1/2, subdivision = (4,.1, 4)
unstableFileID = [savePath, 'P2_unstable_', runID]; % filename to save unstable manifold atlas data
stableFileID = [savePath, 'P1_stable_', runID]; % filename to save stable manifold atlas data
collisionFileID = [savePath, 'P2_P1_connections_', runID]; % filename to save collision data
%% ================================================== SET RUN PARAMETERS  ==================================================
% Manifold run parameters
mu = 0.5; % equal masses means m1 = m2 = 1/2
C = 3; % L4 energy for equal masses
minTau = -5; % stable manifold integration time
maxTau = 5; % unstable manfiold integration time
initialStableArc = [pi/2, pi/2]; % [midPoint, arcRadius] for initial data on P1 collision manifold. [pi/2, pi/2] parameterizes only a half circle
%   to take advantage of antipodal symmetry for equal masses
initialUnstableArc = [pi/2, pi/2]; % duplicate initial data for forward time collision manifold

% Integrator parameters
N = 20; % spatial truncation
M = 40; % temporal truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
tauGuess = 0.1; % initial timestep to scale each Taylor step
%%  ================================================== GROW OR LOAD UNSTABLE MANIFOLD FOR P2  ==================================================
if isfile(unstableFileID) % check if unstable manifold is computed already
    disp('Using existing unstable manifold data')
    
else % compute and save unstable manifold data
    disp('Computing new unstable manifold data')
    initialData = setreginitialdata(mu, 2, N, initialUnstableArc, false); % set up some initial data on a F1 circle
    
    % Integrate the collision manifold for the large primary backwards in time.
    unstableBd = RegCRTBPChart(initialData, 'Taylor', 0, truncation, [mu, C], 2, 'InitialScaling', tauGuess, 'boundary', true);
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasP2F = Atlas(unstableBd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasP2F.LeafStack)
        fprintf('%0.4f \n', min([atlasP2F.LeafStack.TimeSpan]))
        atlasP2F.growboundary('RegTime', @regtime)
    end
    save(stableFileID, 'unstableBd', 'atlasP2F') % save manifold data
end
%%  ================================================== GROW OR LOAD STABLE MANIFOLD FOR P1  ==================================================
if isfile(stableFileID) % check if unstable manifold is computed already
    disp('Using existing stable manifold data')
    
else % compute and save stable manifold data
    disp('Computing new stable manifold data')
    initialData = setreginitialdata(mu, 1, N, initialStableArc, false); % set up some initial data on a F1 circle
    
    % Integrate the collision manifold for the large primary backwards in time.
    stableBd = RegCRTBPChart(initialData, 'Taylor', 0, truncation, [mu, C], 1, 'InitialScaling', -tauGuess, 'boundary', true);
    bdCheck = @(obj, boundaryChart, minTau)boundarycheck(obj, boundaryChart, minTau, 'SubDivideParameter', subDivideParameter);
    atlasP1B = Atlas(stableBd, -tauGuess, bdCheck, @advectioncheck, 'MaxTau', minTau);
    while ~isempty(atlasP1B.LeafStack)
        fprintf('%0.4f \n', min([atlasP1B.LeafStack.TimeSpan]))
        atlasP1B.growboundary('RegTime', @regtime)
    end
    save(stableFileID, 'stableBd', 'atlasP1B') % save manifold data
end
%% ================================================== MINE FOR INTERSECTIONS ==================================================
if isfile(collisionFileID) % check if collisions have already been computed
    disp('Using existing collision data')
else
%     % load stable/unstable manifolds
%     load(stableFileID)
%     atlasB = atlasP1B; % set (stable) backward time atlas variable
%     load(unstableFileID)
%     atlasF = atlasP2F; % set (unstable) forward time atlas variable
    
    % mine for connections by leapfrogging fundamental domains in each atlas
    t0 = 0; % local coordinate to fix in Chart 2 (forward chart) for Newton iteration
    sols = {}; % initialize cell array for connections
    genIdxB = [atlasB.Chart.Generation]; % indices for each fundamental domain
    genIdxF = [atlasF.Chart.Generation]; % indices for each fundamental domain
    genB = -1; % stable atlas generation index
    genF = 0;  % unstable atlas generation index
    genFCurrent = atlasF.generation(genF); % initialize unstable fundamental domain as a vector of charts
    nGenF = length(genFCurrent); % number of charts parameterizing the initial unstable fundamental domain
    
    while genB < atlasB.LastGeneration || genF < atlasF.LastGeneration
        
        % leapfrog to the next stable fundamental domain
        genB = min(genB + 1, atlasB.LastGeneration);
        genBCurrent = atlasB.generation(genB); % parameterize stable fundamental domain as a vector of charts
        fprintf('Mining generations: %d-%d \n',[genB,genF])
        nGenB = length(genBCurrent); % number of charts parameterizing the current stable fundamental domain
        
        % double loop over charts in both current fundamental domains
        for iGenB = 1:nGenB
            for jGenF = 1:nGenF
                chart1 = genBCurrent(iGenB);
                chart2 = genFCurrent(jGenF);
                if isequal(chart1.RegType, chart2.RegType) % check if charts are in the same coordinates
                    ijSols = check4intersection(chart1, chart2, 5);
                    if ~isempty(ijSols) % a connection was found. Append it to the array
                        sols{end+1} = {genBCurrent(iGenB); genFCurrent(jGenF); ijSols};
                    end
                end
            end
        end
        
        % leapfrog to the next unstable fundamental domain
        genF = min(genF + 1, atlasF.LastGeneration);
        genFCurrent = atlasF.generation(genF); % parameterize unstable fundamental domain as a vector of charts
        fprintf('Mining generations: %d-%d \n',[genB,genF])
        nGenF = length(genFCurrent); % number of charts parameterizing the current unstable fundamental domain
        
        % double loop over charts in both current fundamental domains
        for iGenB = 1:nGenB
            for jGenF = 1:nGenF
                chart1 = genBCurrent(iGenB);
                chart2 = genFCurrent(jGenF); 
                if isequal(chart1.RegType, chart2.RegType)
                    ijSols = check4intersection(chart1,chart2,5);
                    if ~isempty(ijSols)
                        sols{end+1} = {genBCurrent(iGenB); genFCurrent(jGenF); ijSols};
                    end
                end
            end
        end
    end
    save(collisionFileID, 'sols') % save collision data
end