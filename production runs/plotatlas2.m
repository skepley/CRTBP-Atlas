%PLOTATLAS2 - Compute plots of the (un)stable atlases computed in final2.m
%   Plot (un)stable manifolds for the following example:
%       Parameters: equal mass
%       Energy: L4
%       Stable manifold: L4
%       Unstable manifold: L4
%       Time: +/- 5 units
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 06-Feb-2020; Last revision: 06-Feb-2020


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

%% ================================================== SET RUN PARAMETERS AND LOAD MANIFOLDS  ==================================================
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

% load stable/unstable manifolds
load l4_equalmass_local_manifolds % get local manifold data
load(stableFileID)
atlasB = atlasL4B; % set (stable) backward time atlas variable
load(unstableFileID)
atlasF = atlasL4F; % set (unstable) forward time atlas variable


%% ================================================== PLOT THE STABLE MANIFOLD ==================================================
close all
figure
hold on


% TEST BED FOR THE CRTBPpatch method

% INPUTS: obj, gridArray, chartRegType, plotRegType plotIdx
obj = atlasF;
globalSpace = linspace(-1, 1, 100);
globalTime = linspace(0, maxTau, 100);
nodeArray = {globalSpace, globalTime};


plotIdx = [1, 3];
f0Color = [0, 1, 0];
f1Color = [1, 0, 0];
f2Color = [0, 0, 1];


CRTBPpatch(atlasF, nodeArray, 0, 0, plotIdx, 'FaceColor', f0Color);
CRTBPpatch(atlasF, nodeArray, 1, 0, plotIdx, 'FaceColor', f1Color);
CRTBPpatch(atlasF, nodeArray, 2, 0, plotIdx, 'FaceColor', f2Color);

%%  ================================================== Function block begins here ==================================================
for iChart = obj.Chart% loop through charts of correct regtype  
    iChart.plotboundary([40,40], 0, plotIdx, 'PlotOptions', {'LineWidth', 1, 'Color', [0, 0, 1]})
end

%  plot primaries and L4
primarySize = 500;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','k')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')
plot3(L4(1), L4(3), L4(2), 'r*')
plot(L4(1), L4(3), 'r*')

axis tight
ax = axis;
dealfig()

