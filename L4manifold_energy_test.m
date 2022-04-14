%L4MANIFOLD_ENERGY_TEST - debug the weird energy drift which appears in the
%integrator
%
%   Other m-files required: none
%   MAT-files required: l4_equalmass_local_manifolds

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 20-Sep-2019; Last revision: 20-Sep-2019

% SET PATH AND SAVE VARIABLES
clear all
% close all
clc
computerPath = pcpath('mac');

% addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath('/Users/shane/Dropbox/Regularisation3bp/CRTBP Atlas'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

%% ================================================== GROW A TEST STABLE MANIFOLD FOR L4 ==================================================
load l4_equalmass_local_manifolds % get local manifold data

% % Run 1: t = 3 with default subdivision strategy
% subDivideParameter = [4, .1, 4]; 
% saveFile = [savePath, 'L4_t3_stable_energy_test_default.mat'];
% maxTau = -3;


% Run 3: t = 3 with strict subdivision strategy
% saveFile = [savePath, 'L4_t3_stable_energy_test_strict.mat'];
% subDivideParameter = [4, .001, 4]; 
% maxTau = -3;


% Run 4: t = 3.25 with default subdivision strategy
% subDivideParameter = [4, .1, 4]; 
% saveFile = [savePath, 'L4_t325_stable_energy_test_default.mat'];
% maxTau = -3.25;

% % Run 5: t = 3.25 with moderate subdivision strategy
% subDivideParameter = [4, .01, 4]; 
% saveFile = [savePath, 'L4_t325_stable_energy_test_moderate.mat'];
% maxTau = -3.25;

% Run 6: t = 3.25 with strict subdivision strategy
% subDivideParameter = [4, .001, 4]; 
% saveFile = [savePath, 'L4_t325_stable_energy_test_strict.mat'];
% maxTau = -3.25;

% Run 7: t = 4 with default subdivision strategy
% subDivideParameter = [4, .1, 4]; 
% saveFile = [savePath, 'L4_t4_stable_energy_test_default.mat'];
% maxTau = -4;

% ====================== JUST COMPUTE A SECTOR ======================
% Run 8: t = 3.5 sector with default subdivision strategy
subDivideParameter = [4, .1, 4]; 
saveFile = [savePath, 'L4_35_stable_energy_test_sector_default.mat'];
maxTau = -3.5;

% % Run 9: t = 3.5 sector with strict subdivision strategy
% subDivideParameter = [4, .001, 4]; 
% saveFile = [savePath, 'L4_35_stable_energy_test_sector_strict.mat'];
% maxTau = -3.5;

% % Run 9: t = 3.5 sector with ultra-strict subdivision strategy
% subDivideParameter = [4, .0001, 4]; 
% saveFile = [savePath, 'L4_35_stable_energy_test_sector_ultrastrict.mat'];
% maxTau = -3.5;
% % tau = -.1;
% 
% % % Run 10: t = 3.5 sector with strict subdivision strategy and test step =
% % 0.01;
% subDivideParameter = [4, .001, 4]; 
% saveFile = [savePath, 'L4_35_newtau_stable_energy_test_sector_strict.mat'];
% maxTau = -3.5;
% tau = -.01;


% Integrator and collision manifold parameters
basis = 'Taylor';
initialTime = 0; % always 0 for autonomous system
isValid = false; % validate computations
N = 30; % spatial truncation
truncation = [20, N]; % time and space truncation spaces
regType = 0;
parameter = mu;

% load existing manifold data or compute and save new data
if isfile(saveFile)
    load(saveFile)
else  % compute the stable manifold
    
    % lift local stable parameterization into phase space
   
    % pick off coordinates
    XsLocal = mid(As);
    PsLocal = mid(Bs);
    YsLocal = mid(Cs);
    QsLocal = mid(Ds);
    RsLocal = mid(Rs);
    SsLocal = mid(Ss);
    
    numLineSeg = 8;
    
    % take equally spaced nodes on the circular sector

    linearSpan = [-0.3, -0.1]; % Take only the spatialSpan interval
    sectorSpan = pi*linearSpan + pi; % corresponding sector
    theta = linspace(sectorSpan(1), sectorSpan(2), numLineSeg + 1);
    nodes = exp(1i*theta);
    globalMapCoef = [sectorSpan', [1;1]]\[-1;1]; % [m,b] where T(x) = mx + b is the linear bijection from [a,b] to [-1,1]
    globalSpatialSpan = polyval(globalMapCoef, theta);
    
    % GLOBALSPATIALSPAN should just be given by linspace(-1,1,numLinSeg).
    % The rest of this code is correct but needlessly complicated. 


    % Take equally spaced nodes on the entire unit circle
    %     theta = linspace(0, 2*pi, numLineSeg + 1);
    %     nodes = exp(1i*theta);
    %     globalSpatialSpan = (theta - pi)/pi;

    
    % lift parameterization
    clear stableBd;
    for k = 1:numLineSeg
        
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
        localStableBd = real([X0(1:truncation(2));P0(1:truncation(2));Y0(1:truncation(2));Q0(1:truncation(2));R0(1:truncation(2));S0(1:truncation(2))]);
        kSpatialSpan = globalSpatialSpan(k:k+1); % get coordinates in interval [-1,1]
        stableBd(k) =  RegCRTBPChart(localStableBd, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
        stableBd(k).SpatialSpan = kSpatialSpan;
    end
    
    % Integrate the stable manifold backwards in time.
    bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
    atlasL4B = Atlas(stableBd, tau, bdCheck, @advectioncheck, 'MaxTau', maxTau);
    while ~isempty(atlasL4B.LeafStack)
        fprintf('%0.4f \n', min([atlasL4B.LeafStack.TimeSpan]))
        atlasL4B.growboundary('RegTime', @regtime)
    end
    save(saveFile)
end
check_atlas_energy




