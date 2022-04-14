%UNFUCK_L1_ALGEBRA - Fix ell^1 algebra on initial conditions which appears
% wrong somehow. It is off by around 10^-1 (see Test 1). But it does not always cause bad orbits (see Test 2).
%
%
%   Other m-files required: none
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 10-Jan-2021; Last revision: 10-Jan-2021

clear all
clc
close all




%%  COMPUTE CONTINUATION OF A SINGLE INITIAL CONDITION
% % CHECK F0 AGAINST RUNGE KUTTA
% parameter = U.Parameter;
% maxTau = U.TimeSpan(2);
% truncation = U.Truncation;
% basis = 'Taylor';
% initialTime = 0;
% regType = 0;
% tauGuess = 0.1;
% subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
% hotSwap = false; % specify whether charts should swap based on the ideal domain map or not
% bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter,...
%     'HotSwap', hotSwap);
% 
% % set up some initial data
% initialData = [X0.Coefficient; P0.Coefficient; Y0.Coefficient; Q0.Coefficient];
% bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
% 
% % Integrate multiple timesteps
% A0 = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
% while ~isempty(A0.LeafStack)
%     A0.growboundary()



%% ================================================== TEST 4 ==================================================
% lets look at the ell^1 algebra on the initial data

% COPIED FROM FINAL1.M
% set paths for accessing and saving data
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

% ================================================== SET RUN PARAMETERS  ==================================================
% Manifold run parameters
mu = 0.5; % equal masses means m1 = m2 = 1/2
C = 3; % L4 energy for equal masses
maxTau = 5; % unstable manfiold integration time
minTau = -5;
numUnstableSegment = 8; % number of initial subdivisions for unstable boundary
regTime = @taylorregtime;  % set method of integration for regularizing time

% Integrator parameters
N = 20; % spatial truncation
M = 40; % temporal truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
tauGuess = 0.1; % initial timestep to scale each Taylor step

% % local stable/unstable preimage mappings (curry run parameters)
% unstableLocalMap = @(s)localunstablecoordinates(numUnstableSegment, s);
% stableLocalMap = @(s)localstablecoordinates(initialStableArc, s);

% first lift local L4 unstable manifold into phase space and get boundary
load newl4_equalmass_local_manifolds % get local manifold data and pick off coordinates

% check unstable manifold
XuLocal = mid(Au);
PuLocal = mid(Bu);
YuLocal = mid(Cu);
QuLocal = mid(Du);
RuLocal = mid(Ru);
SuLocal = mid(Su);

% % check stable manifold
% XuLocal = mid(As);
% PuLocal = mid(Bs);
% YuLocal = mid(Cs);
% QuLocal = mid(Ds);
% RuLocal = mid(Rs);
% SuLocal = mid(Ss);

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

% plot orbits fron initial points on manifold boundary
parameter = mu;
% close all
t = linspace(0, maxTau, 100); % plotting time points
s = linspace(-1, 1,10); % some spatial sample points
regType = 0;
VF = @(t,x)rk45regvectorfield(t,x,parameter, regType);
tSpan = [0, minTau];
figure;
hold on
for k = 1:length(unstableBd)
    init_U = unstableBd(k); % this is the initial data for U used in the computation

    for j = 1:numel(s)
        % plot rk45 orbit using after evaluating initial parameterization
        sj = s(j);
        ob = cell2mat(init_U.eval(sj));
        plot_orbit(VF, tSpan, ob', [1,3], 'PlotOptions', {'g'})
        
        % plot rk45 orbit using only first four coordinates
        u0 = ob(1:4)';
%         r0 = 1./sqrt((u0(1) - mu).^2 + u0(3).^2);
%         s0 = 1./sqrt((u0(1) + 1 - mu).^2 + u0(3).^2);
%         u0 = [u0; r0; s0];
%         norm(u0 - ob)
        plot_orbit(VF, tSpan, addAutoDiff(u0, mu), [1,3], 'PlotOptions', {'r-.'})
    end
end

%% ================================================== TEST 5 ==================================================
% Check if the precision loss is due to the local parameterization
% take equally spaced nodes on the circle


% theta = linspace(0, 2*pi, 100);
% z = exp(1i*theta).';
% zBar = conj(z);
% xLoc = Scalar(mid(Au), 'Taylor');
% pLoc = Scalar(mid(Bu), 'Taylor');
% yLoc = Scalar(mid(Cu), 'Taylor');
% qLoc = Scalar(mid(Du), 'Taylor');
% rLoc = Scalar(mid(Ru), 'Taylor');
% sLoc = Scalar(mid(Su), 'Taylor');
% Uloc = Chart();
% Uloc.Coordinate = [xLoc; pLoc; yLoc; qLoc; rLoc; sLoc];
% Uloc.Dimension = [2, 6];
% bdLocCellEval = Uloc.eval([z, zBar]);
% bdLocEval = [bdLocCellEval{:}];
% 
% figure
% hold on
% plot(bdLocEval(:,1), bdLocEval(:,3))
% plot(L4(1), L4(3), 'r*') % plot to make sure this is the boundary of L4 local manifold
% rChk = 1./sqrt((bdLocEval(:, 1) - mu).^2 + bdLocEval(:, 3).^2);
% rChk2 = sqrt((bdLocEval(:, 1) - mu).^2 + bdLocEval(:, 3).^2);
% sChk = 1./sqrt((bdLocEval(:, 1) + 1 - mu).^2 + bdLocEval(:, 3).^2);
% sChk2 = sqrt((bdLocEval(:, 1) + 1 - mu).^2 + bdLocEval(:, 3).^2);
% 
% mean(abs(rChk - bdLocEval(:, 5)))
% mean(abs(sChk - bdLocEval(:, 6)))
% mean(abs(rChk2 - bdLocEval(:, 5)))
% mean(abs(sChk2 - bdLocEval(:, 6)))

function uFull = addAutoDiff(u, mu)
r0 = 1./sqrt((u(1) - mu).^2 + u(3).^2);
s0 = 1./sqrt((u(1) + mu).^2 + u(3).^2);
uFull = [u; r0; s0];
end


