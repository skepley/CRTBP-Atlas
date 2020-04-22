%FINAL2 - Production run for homoclinic connections to L4 at equal masses.
%   Connecting orbits for the following example:
%       Parameters: equal mass
%       Energy: L4
%       Forward target: L4
%       Backward target: L4
%       Time: +/- 5 units
%
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 31-Dec-2019; Last revision: 31-Dec-2019

%% ================================================== SET DATA PATHS AND FILENAMES  ==================================================
clear all
close all
clc

% set paths for accessing and saving data
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/SeqDE']))
addpath(genpath('/Users/shane/Dropbox/Regularisation3bp/CRTBP Atlas'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

% set file names to save manifold and collision data
runID = 't5_CL4_4020_equalmass_default.mat'; % t = 5, C = 3, (M,N) = (40,20), mu = 1/2, subdivision = (4,.1, 4)
unstableFileID = [savePath, 'L4_unstable_', runID]; % filename to save unstable manifold atlas data
stableFileID = [savePath,'L4_stable_', runID]; % filename to save stable manifold atlas data
collisionFileID = [savePath, 'L4_L4_connections_', runID]; % filename to save collision data
%% ================================================== SET RUN PARAMETERS  ==================================================
% Manifold run parameters
mu = 0.5; % equal masses means m1 = m2 = 1/2
C = 3; % L4 energy for equal masses
minTau = -5; % stable manifold integration time
maxTau = 5; % unstable manfiold integration time
numUnstableSegment = 8; % number of initial subdivisions for unstable boundary
numStableSegment = 8; % number of initial subdivisions for stable boundary

% Integrator parameters
N = 20; % spatial truncation
M = 40; % temporal truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
tauGuess = 0.1; % initial timestep to scale each Taylor step
%%  ================================================== GROW OR LOAD UNSTABLE MANIFOLD FOR L4  ==================================================
if isfile(unstableFileID) % check if unstable manifold is computed already
    disp('Using existing unstable manifold data')
    
else % compute and save unstable manifold data
    disp('Computing new unstable manifold data')

    % first lift local L4 unstable manifold into phase space and get boundary
    load l4_equalmass_local_manifolds % get local manifold data and pick off coordinates
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
    
    % Integrate the unstable manifold forward in time.
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasL4F = Atlas(unstableBd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasL4F.LeafStack)
        fprintf('%0.4f \n', min([atlasL4F.LeafStack.TimeSpan]))
        atlasL4F.growboundary('RegTime', @regtime)
    end
    save(unstableFileID, 'unstableBd', 'atlasL4F') % save manifold data
end
%%  ================================================== GROW OR LOAD STABLE MANIFOLD FOR L4  ==================================================
if isfile(stableFileID) % check if stable manifold is computed already
    disp('Using existing stable manifold data')
    
else % compute and save stable manifold data
    disp('Computing new stable manifold data')

    % first lift local L4 stable manifold into phase space and get boundary
    load l4_equalmass_local_manifolds % get local manifold data and pick off coordinates
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
    clear stableBd;
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
        stableBd(k) =  RegCRTBPChart(localStableBd, 'Taylor', 0, truncation, mu, 0, 'InitialScaling', -tauGuess, 'boundary', true);
        stableBd(k).SpatialSpan = kSpatialSpan;
    end
    
    % Integrate the stable manifold backward in time.
    bdCheck = @(obj, boundaryChart, minTau)boundarycheck(obj, boundaryChart, minTau, 'SubDivideParameter', subDivideParameter);
    atlasL4B = Atlas(stableBd, -tauGuess, bdCheck, @advectioncheck, 'MaxTau', minTau);
    while ~isempty(atlasL4B.LeafStack)
        fprintf('%0.4f \n', min([atlasL4B.LeafStack.TimeSpan]))
        atlasL4B.growboundary('RegTime', @regtime)
    end
    save(stableFileID, 'stableBd', 'atlasL4B') % save manifold data
end

%% ================================================== MINE FOR INTERSECTIONS ==================================================
if isfile(collisionFileID) % check if collisions have already been computed 
    disp('Using existing collision data')
else
    % load stable/unstable manifolds
    load(stableFileID)
    atlasB = atlasL4B; % set (stable) backward time atlas variable
    load(unstableFileID)
    atlasF = atlasL4F; % set (unstable) forward time atlas variable
    
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