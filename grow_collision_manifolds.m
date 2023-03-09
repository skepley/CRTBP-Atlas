%GROW_COLLISION_MANIFOLDS - Look for a homoclinic connection for a collision manifold
%
%   Description:
%       GROW_COLLISION_MANIFOLDS description
%
%   Output:
%       GROW_COLLISION_MANIFOLDS output
%
%   Other m-files required: none
%   MAT-files required: none
% 
%   See also: GROW_L4_MANIFOLDS, MINE_CONNECTIONS

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 13-Jun-2019; Last revision: 2-Oct-2019 

clear all
close all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath('/Users/shane/Dropbox/Regularisation3bp/CRTBP Atlas'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

%% ================================================== SET UP PARAMETERS FOR DIFFERENT RUNS ==================================================
% Integrator and collision manifold parameters

% RUN 1: L4 energy in the equal mass system 
% mu = 1/2; % small mass primary
% M = 40; % temporal truncation
% N = 20; % spatial truncation 
% C = 3; % L4 energy 
% maxTime = 5;
% subDivideParameter = [4, .1, 4]; % default is (4, .1, 4)
% initialDataArc = [pi/2, pi/2]; % [midPoint, arcRadius] for arc of initial data. Set to [pi/2, pi/2] to take advantage of antipodal symmetry for symmetric masses

% runID = 't5_CL4_4020_eqmass_default.mat'; % t = 5, C = 3, (M,N) = (40,20), mu = 1/2, subdivision = (4,.1, 4)



% RUN 2: L4 energy in the earth-moon system  
mu = 1/81; % small mass primary
M = 40; % temporal truncation 
N = 20; % spatial truncation 
C = 3.6725805575177963; % L4 energy 
maxTime = 5;
subDivideParameter = [4, .1, 4]; % default is (4, .1, 4)
initialDataArc = [pi/2, pi]; % [midPoint, arcRadius] for arc of initial data. Set to [pi/2, pi] for nonsymmetric cases 

runID = 't5_CL4_4020_earthmoon_default.mat'; % t = 5, C = 3, (M,N) = (40,20), mu = 1/2, subdivision = (4,.1, 4)



% set the remaining parameters
basis = 'Taylor';
truncation = [M, N]; % time and space truncation spaces
initialTime = 0; % always 0 for autonomous system
isValid = false; % validate computations
parameter = [mu,C];


%% ================================================== GROW P1 MANIFOLDS ==================================================
saveFile = [savePath,'P1_stable_',runID]; % filename to save data
if isfile(saveFile) % avoid overwriting saved manifold data
    disp('P1 manifolds with this filename already exist')
    
else % Integrate the collision manifold for the large primary backwards in time. 
    
    tic
    % set initial data
    regType = 1;
    initialDataArc = [pi/2, pi/2]; % [midPoint, arcRadius] for arc of initial data. Set to [pi/2, pi/2] to take advantage of antipodal symmetry for symmetric masses
    initialData = setreginitialdata(mu, regType, N, initialDataArc, false); % set up some initial data on a F1 circle
    
    % Integrate backward
    tau = -.1;
    maxTau = -maxTime;
    bdP1B = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasP1B = Atlas(bdP1B, tau, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasP1B.LeafStack)
        fprintf('%0.4f \n', max([atlasP1B.LeafStack.TimeSpan]))
        atlasP1B.growboundary('RegTime', @regtime)
    end
    save(saveFile)
    clear atlasP1B; clear bdP1B;
    
    % ================================================== GROW P1 UNSTABLE MANIFOLD ==================================================
    saveFile = [savePath,'P1_unstable_',runID]; % filename to save data
    
    % Integrate the collision manifold for the large primary forwards in time.
    tau = .1;
    maxTau = maxTime;
    bdP1F = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasP1F = Atlas(bdP1F, tau, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasP1F.LeafStack)
        fprintf('%0.4f \n', min([atlasP1F.LeafStack.TimeSpan]))
        atlasP1F.growboundary('RegTime', @regtime)
    end
    save(saveFile)
    clear atlasP1F; clear bdP1F;
    toc
end


%% ================================================== GROW P2 MANIFOLDS ==================================================
saveFile = [savePath,'P2_stable_',runID]; % filename to save data
if isfile(saveFile) % avoid overwriting saved manifold data
    disp('P2 manifolds with this filename already exist')
    
else % Integrate the collision manifold for the large primary backwards in time.
    tic
    
    % set initial data
    regType = 2;
    initialData = setreginitialdata(mu, regType, N, initialDataArc, false); % set up some initial data on a F2 circle
    
    % Integrate the collision manifold for the SMALL primary backwards in time.
    tau = -.1;
    maxTau = -maxTime;
    bdP2B = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasP2B = Atlas(bdP2B, tau, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasP2B.LeafStack)
        fprintf('%0.4f \n', max([atlasP2B.LeafStack.TimeSpan]))
        atlasP2B.growboundary('RegTime', @regtime)
    end
    save(saveFile)
    clear atlasP2B; clear bdP2B;
    
    % ================================================== GROW P2 UNSTABLE MANIFOLD ==================================================
    saveFile = [savePath,'P2_unstable_',runID]; % filename to save data
    
    % Integrate the collision manifold for the SMALL primary forwards in time.
    tau = .1;
    maxTau = maxTime;
    bdP2F = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasP2F = Atlas(bdP2F, tau, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasP2F.LeafStack)
        fprintf('%0.4f \n', min([atlasP2F.LeafStack.TimeSpan]))
        atlasP2F.growboundary('RegTime', @regtime)
    end
    save(saveFile)
    clear atlasP2F; clear bdP2F;
    toc
end

%% ================================================== % MINE FOR INTERSECTIONS  ==================================================

% THIS CODE HAS BEEN MOVED TO THE SCRIPT NAMED MINE_CONNECTIONS.M

return
%% TESTING SURF PLOTS
close all

% plot the manifold in F0 coordinates
tFinal = atlasB.MaxTau;
nPoint = [3,3];
s = linspace(-1,1,nPoint(1));
t = linspace(0, 1, nPoint(2));
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];

figure
hold on
for jChart = atlasB.Chart
    % map evaluations to F0 coordinates
    if isequal(jChart.RegType, 1)
        [X1,P1,Y1,Q1,~] = jChart.eval(evalData);
        F0Data = CRTBP2reg([X1,P1,Y1,Q1],mu,-1);
    elseif isequal(jChart.RegType, 2)
        [X2,P2,Y2,Q2,~] = jChart.eval(evalData);
        F0Data = CRTBP2reg([X2,P2,Y2,Q2],mu,-2);
    else
        [X0, P0, Y0, Q0, ~, ~] = jChart.eval(evalData);
        F0Data = [X0, P0, Y0, Q0];
    end
    
    X = reshape(F0Data(:,1), length(s), []);
    Y = reshape(F0Data(:,3), length(s), []);
    Z = reshape(F0Data(:,2), length(s), []);
    surf(X,Y,Z, 'FaceColor', [0,1,0])
end

tFinal = atlasF.MaxTau;
nPoint = [3,3];
s = linspace(-1,1,nPoint(1));
t = linspace(0, 1, nPoint(2));
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];

% figure
% hold on
for jChart = atlasF.Chart
    % map evaluations to F0 coordinates
    if isequal(jChart.RegType, 1)
        [X1,P1,Y1,Q1,~] = jChart.eval(evalData);
        F0Data = CRTBP2reg([X1,P1,Y1,Q1],mu,-1);
    elseif isequal(jChart.RegType, 2)
        [X2,P2,Y2,Q2,~] = jChart.eval(evalData);
        F0Data = CRTBP2reg([X2,P2,Y2,Q2],mu,-2);
    else
        [X0, P0, Y0, Q0, ~, ~] = jChart.eval(evalData);
        F0Data = [X0, P0, Y0, Q0];
    end
    
    X = reshape(F0Data(:,1), length(s), []);
    Y = reshape(F0Data(:,3), length(s), []);
    Z = reshape(F0Data(:,2), length(s), []);
    surf(X,Y,Z, 'FaceColor', [1,0,0])
end

axis([mu-1-.5, mu+.5, -1, 1, -3, 3])
view(-121,16)
dealfig()



