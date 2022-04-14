%MINE_CONNECTIONS - Mine for connections between a forward and a backward time manifold 

%
%   Description:
%       MINE_CONNECTIONS description
%
%   Output:
%       MINE_CONNECTIONS output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 07-Jul-2019; Last revision: 2-Oct-2019
 

clear all
close all 
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath('/Users/shane/Dropbox/Regularisation3bp/CRTBP Atlas'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

%% ================================================== % LOAD STABLE/UNSTABLE MANIFOLDS ==================================================
% runID = 't5_CL4_4020_eqmass_default.mat'; % t = 5, C = L4 energy, (M,N) = (40,20), mu = 1/2, subdivision = (4,.1, 4)
runID = 't5_CL4_4020_earthmoon_default.mat'; % t = 5, C = L4 energy, (M,N) = (40,20), mu = 1/81, subdivision = (4,.1, 4)

% ========================= LOAD A BACKWARD TIME MANIFOLD =========================
stableAtlasID = 'P1_stable_';
load([savePath,stableAtlasID,runID])
% set stable atlas variable
atlasB = atlasP1B; 
% atlasB = atlasP2B;

% =========================LOAD A FORWARD TIME MANIFOLD =========================
unstableAtlasID = 'P2_unstable_';
load([savePath,unstableAtlasID,runID])
% set unstable atlas variable
% atlasF = atlasP1F;
atlasF = atlasP2F;


%% ================================================== MINE FOR INTERSECTIONS ==================================================
saveCollisionFile = [unstableAtlasID(1:3), stableAtlasID(1:2), 'connections_', runID]; % set filename to save collision data

t0 = 0; % local coordinate to fix in Chart 2 for Newton iteration
sols = {};
genIdxB = [atlasB.Chart.Generation];
genIdxF = [atlasF.Chart.Generation];
genB = -1; % bakward collision atlas generation index
genF = 0;  % forward collision atlas generation index
genFCurrent = atlasF.generation(genF);
nGenF = length(genFCurrent);

while genB < atlasB.LastGeneration || genF < atlasF.LastGeneration
    genB = min(genB + 1, atlasB.LastGeneration);
    genBCurrent = atlasB.generation(genB);
    disp([genB,genF])
    nGenB = length(genBCurrent);
    
    % loop over current atlas1 generation
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
    
    % leapfrog to a fundamental domain for the next generation of atlas 2
    genF = min(genF + 1, atlasF.LastGeneration);
    genFCurrent = atlasF.generation(genF);
    disp([genB,genF])
    nGenF = length(genFCurrent);
    
    % loop over backward gen
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

save([savePath, saveCollisionFile], 'sols', 'C', 'mu') % save collision data


return

%% ================================================== CODE FROM BEFORE COMPUTER WIPE ==================================================
% OLD MANIFOLD SAVE DATA
% load([savePath, 'P1_equalmass_L4energy_t5backward'])
% atlasB = atlasL4B;


% OLD MANIFOLD SAVE DATA
% load([savePath, 'P1_equalmass_L4energy_t5forward'])
% atlasF = atlasL4F;

% load([savePath, 'P2_equalmass_L4energy_t5forward'])
% atlasF = atlasP2F;

% load([savePath, 'P2_earthmoon_L4energy_t5forward'])
% atlasF = atlasP2F;



% save([savePath, 'P1P2_collision_heteroclinics_equalmass_L4energy'], 'sols', 'C', 'mu')
% save([savePath, 'P1P2_collision_heteroclinics_earthmoon_L4energy'], 'sols', 'C', 'mu')
% save([savePath, 'P1_collision_homoclinics_equalmass_L4energy'], 'sols', 'C', 'mu')
% save([savePath, 'L4_homoclinics_equalmass'], 'sols', 'C', 'mu')



