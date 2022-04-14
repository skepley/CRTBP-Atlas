%FINAL2 - Production run for homoclinic connections to L4 at equal masses.
%   Connecting orbits for the following example:
%       Parameters: equal mass
%       Energy: L4
%       Forward target: L4
%       Backward target: L4
%       Time: +/- 7 units
%
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 31-Dec-2019; Last revision: 10-Mar-2021

%% ================================================== SET DATA PATHS AND FILENAMES  ==================================================
clear all
close all
clc 


% addpath('/users/shane/dropbox/matlab/') % add this for the path on Elena's computer 
% addpath(genpath('/users/shane/downloads/intlab_V12'))

% set paths for accessing and saving data
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

% set file names to save manifold and collision data
% runID = 't7_CL4_4021_equalmass_default'; % t = 7, C = 3, (M,N) = (40,21), mu = 1/2, subdivision = (4,.1, 4)
% runID = 't5_CL4_4021_equalmass_default'; % t = 5, C = 3, (M,N) = (40,21), mu = 1/2, subdivision = (4,.1, 4)
runID = 't7_CL4_3111_equalmass_default'; % t = 7, C = 3, (M,N) = (40,21), mu = 1/2, subdivision = (4,.1, 4)
unstableFileID = ['L4_unstable_', runID]; % filename to save unstable manifold atlas data
stableFileID = ['L4_stable_', runID]; % filename to save stable manifold atlas data
collisionFileID = [savePath, 'L4_L4_connections_', runID, '.mat']; % filename to save collision data


%% ================================================== SET RUN PARAMETERS  ==================================================
% Manifold run parameters
mu = 0.5; % equal masses means m1 = m2 = 1/2
C = 3; % L4 energy for equal masses
minTau = -7; % stable manifold integration time
maxTau = 7; % unstable manfiold integration time
numUnstableSegment = 8; % number of initial subdivisions for unstable boundary
numStableSegment = 8; % number of initial subdivisions for stable boundary
regTime = @taylorregtime;  % set method of integration for regularizing time

% Integrator parameters
N = 11; % spatial truncation
M = 31; % temporal truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
tauGuess = 0.1; % initial timestep to scale each Taylor step

% local stable/unstable preimage mappings (curry run parameters)
unstableLocalMap = @(s)localunstablecoordinates(numUnstableSegment, s);
stableLocalMap = @(s)localstablecoordinates(numStableSegment, s);

    
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



%%  ================================================== GROW OR LOAD STABLE MANIFOLD FOR L4  ==================================================
if exist([savePath, stableFileID]) || exist([savePath, stableFileID, '.mat']) % check if unstable manifold is computed already
    disp('Using existing stable manifold data')
    atlasL4B = chunkyload(savePath, stableFileID);
    
else % compute and save stable manifold data
    disp('Computing new stable manifold data')
    stableBd = L4_local_stable_boundary('newL4_equalmass_local_manifolds', numStableSegment, truncation, mu, tauGuess);
    
    % Integrate the stable manifold backward in time.
    bdCheck = @(obj, boundaryChart, minTau)boundarycheck(obj, boundaryChart, minTau, 'SubDivideParameter', subDivideParameter);
    atlasL4B = RegCRTBPAtlas(stableBd, -tauGuess, bdCheck, @advectioncheck, 'MaxTau', minTau);
    while ~isempty(atlasL4B.LeafStack)
        atlasL4B.growboundary('RegTime', regTime)
    end
    chunkysave(atlasL4B, savePath, stableFileID) % save manifold data
end
% fprintf('Stable Atlas Size: %d \n', numel(atlasL4B.Chart))


%% ================================================== MINE FOR INTERSECTIONS ==================================================
if isfile(collisionFileID) % check if collisions have already been computed
    disp('Using existing collision data')
    load(collisionFileID)
    
else
    atlasB = atlasL4B; % set (stable) backward time atlas variable
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
    
    %%
    while genB < atlasB.LastGeneration || genF < atlasF.LastGeneration
        
        % leapfrog to the next stable fundamental domain
        if genB < atlasB.LastGeneration
            %             genB = min(genB + 1, atlasB.LastGeneration);
            genB = genB + 1;
            genBCurrent = atlasB.generation(genB); % parameterize stable fundamental domain as a vector of charts
            fprintf('Mining generations: %d-%d \n',[genB,genF])
            nGenB = length(genBCurrent); % number of charts parameterizing the current stable fundamental domain
            
            % double loop over charts in both current fundamental domains
            for iGenB = 1:nGenB
                for jGenF = 1:nGenF
                    chart1 = genBCurrent(iGenB);  % stable chart
                    chart2 = genFCurrent(jGenF);  % unstable chart
                    [isTrue, ijSols] = check4intersection(chart1, chart2, 5);
                    if ~isempty(ijSols) % a connection was found. Append it to the array
                        sols{end+1} = {chart1; chart2; ijSols; isTrue};
                    end
                end
            end
        end
        
        % leapfrog to the next unstable fundamental domain
        if genF < atlasF.LastGeneration
            genF = genF + 1;
            %         genF = min(genF + 1, atlasF.LastGeneration);
            genFCurrent = atlasF.generation(genF); % parameterize unstable fundamental domain as a vector of charts
            fprintf('Mining generations: %d-%d \n',[genB,genF])
            nGenF = length(genFCurrent); % number of charts parameterizing the current unstable fundamental domain
            
            % double loop over charts in both current fundamental domains
            for iGenB = 1:nGenB
                for jGenF = 1:nGenF
                    chart1 = genBCurrent(iGenB);  % stable chart
                    chart2 = genFCurrent(jGenF);  % unstable chart
                    [isTrue, ijSols] = check4intersection(chart1,chart2,5);
                    if ~isempty(ijSols)
                        sols{end+1} = {chart1; chart2; ijSols; isTrue};
                    end
                end
            end
        end
    end
    save(collisionFileID, 'sols') % save collision data
end


% %% MINING ALL PAIRS FROM EACH ATLAS TO VERIIFY WE HAVE THE CONNECTION
% sols2 = {};
% for i = 1:atlasB.Size
%     fprintf('Stable Chart %d \n', i)
%     for j = 1:atlasF.Size
%         chart1 = atlasB.Chart(i); % stable chart
%         chart2 = atlasF.Chart(j); % unstable chart
%         if isequal(chart1.RegType, chart2.RegType)
%             ijSols = check4intersection(chart1,chart2,5);
%             if ~isempty(ijSols)
%                 sols2{end+1} = {chart1; chart2; ijSols};
%             end
%         end
%     end
% end
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

function localParmCoordinates = localstablecoordinates(nSegment, globalSpace)
% LOCALSTABLECOORDINATES Local manifold coordinate mapping a global spatial coordinate to its preimage under the
% local parameterization of the stable manifold to L4.

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