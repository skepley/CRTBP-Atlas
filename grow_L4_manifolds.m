%GROW_L4_MANIFOLDS - Compute (un)stable manifolds for the equilibrium at L4
%
%   Description:
%       GROW_L4_MANIFOLDS description
%
%   Output:
%       GROW_L4_MANIFOLDS output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: GROW_COLLISION_MANIFOLDS, MINE_CONNECTIONS 
 
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 12-Jun-2019; Last revision: 2-Oct-2019

clear all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath('/Users/shane/Dropbox/Regularisation3bp/CRTBP Atlas'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

%% ================================================== MAKE SOME CHOICES ==================================================
load l4_equalmass_local_manifolds % get local manifold data

% Integrator and collision manifold parameters
basis = 'Taylor';
initialTime = 0; % always 0 for autonomous system
isValid = false; % validate computations
N = 30; % spatial truncation
truncation = [20, N]; % time and space truncation spaces
regType = 0;
parameter = mu;

%% STABLE MANIFOLD
% lift local stable parameterization into phase space
tau = -.1;
maxTau = -3;

% pick off coordinates
XsLocal = mid(As);
PsLocal = mid(Bs);
YsLocal = mid(Cs);
QsLocal = mid(Ds);
RsLocal = mid(Rs);
SsLocal = mid(Ss);

% take equally spaced nodes on the circle
numLineSeg = 8;
theta = linspace(0, 2*pi, numLineSeg + 1);
nodes = exp(1i*theta);
globalSpatialSpan = (theta - pi)/pi;

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


return
% Integrate the stable manifold backwards in time.
subDivideParameter = [4, .001, 4]; % default is [4, .1, 4]
bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
atlasL4B = Atlas(stableBd, tau, bdCheck, @advectioncheck, 'MaxTau', maxTau);
while ~isempty(atlasL4B.LeafStack)
    fprintf('%0.4f \n', min([atlasL4B.LeafStack.TimeSpan]))
    atlasL4B.growboundary('RegTime', @regtime)
end
check_atlas_energy
return
save([savePath, 'L4_equalmass_t5backward'])
clear atlasL4B; clear stableBd;

%% UNSTABLE MANIFOLD
% lift local unstable parameterization into phase space
tau = .1;
maxTau = 5;

% pick off coordinates
XuLocal = mid(Au);
PuLocal = mid(Bu);
YuLocal = mid(Cu);
QuLocal = mid(Du);
RuLocal = mid(Ru);
SuLocal = mid(Su);

% take equally spaced nodes on the circle
numLineSeg = 8;
theta = linspace(0, 2*pi, numLineSeg + 1);
nodes = exp(1i*theta);
globalSpatialSpan = (theta - pi)/pi;

% lift parameterization
clear unstableBd;
for k = 1:numLineSeg
    
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
    localUnstableBd = real([X0(1:truncation(2));P0(1:truncation(2));Y0(1:truncation(2));Q0(1:truncation(2));R0(1:truncation(2));S0(1:truncation(2))]);
    kSpatialSpan = globalSpatialSpan(k:k+1); % get coordinates in interval [-1,1]
    unstableBd(k) =  RegCRTBPChart(localUnstableBd, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
    unstableBd(k).SpatialSpan = kSpatialSpan;
end

% Integrate the stable manifold backwards in time.
subDivideParameter = [3, .01, 4]; % default is [4, .05, 4]
bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
atlasL4F= Atlas(unstableBd, tau,  bdCheck, @advectioncheck, 'MaxTau', maxTau);
while ~isempty(atlasL4F.LeafStack)
    fprintf('%0.4f \n', min([atlasL4F.LeafStack.TimeSpan]))
    atlasL4F.growboundary('RegTime', @regtime)
end
save([savePath, 'L4_equalmass_t5forward'])
clear atlasL4F; clear unstableBd;







return

%% TESTING SURF PLOTS

% plot the manifold in F0 coordinates
s = linspace(-1,1,10);
tFinal = atlasL4F.MaxTau;
t = linspace(0, 1, 10);
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];
close all
figure
hold on

for jChart = atlasL4F.Chart
    
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
    surf(X,Y,Z)
    shading interp
end
colorbar



return
%% PLOT THE MANIFOLDS
close all
figure
hold on

% plot primaries and L4
primarySize = 500;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','b')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')
plot3(L4(1), L4(3), L4(2), 'r*')


if exist('atlasL4F')
    
    % plot the manifold in F0 coordinates
    s = linspace(-1,1,100);
    tFinal = atlasL4F.MaxTau;
    t = linspace(0, tFinal, 500);
    [S,T] = meshgrid(s,t);
    evalData = [reshape(S,[],1), reshape(T,[],1)];
    
    face = [];
    vertex = [];
    colorData = [];
    
    for jChart = atlasL4F.Chart(1:end)
        bdVertices = [jChart.SpatialSpan(1), jChart.TimeSpan(1); jChart.SpatialSpan(1), jChart.TimeSpan(2);... % add 4 corners of domain
            jChart.SpatialSpan(2), jChart.TimeSpan(1); jChart.SpatialSpan(2), jChart.TimeSpan(2)];
        data = [jChart.intersectdomain(evalData);bdVertices];
        
        % map evaluations to F0 coordinates
        if isequal(jChart.RegType, 1)
            [X1,P1,Y1,Q1,~] = jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true);
            F0Data = CRTBP2reg([X1,P1,Y1,Q1],mu,-1);
            jColorData = 0.1*ones(size(data,1),1);
        elseif isequal(jChart.RegType, 2)
            [X2,P2,Y2,Q2,~] = jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true);
            F0Data = CRTBP2reg([X2,P2,Y2,Q2],mu,-2);
            jColorData = 0.5*ones(size(data,1),1);
        else
            [X0, P0, Y0, Q0,~,~] = jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true);
            F0Data = [X0, P0, Y0, Q0];
            jColorData = 0.9*ones(size(data,1),1);
            
        end
        
        X0 = F0Data(:,1);
        P0 = F0Data(:,2);
        Y0 = F0Data(:,3);
        Q0 = F0Data(:,4);
        jVertex = [mid(X0),mid(Y0),mid(P0)];
        jFace = delaunay(data(:,1),data(:,2)); % index triples for delaunay triangulation in space-time.
        vertexCount = size(vertex,1);
        face = [face;jFace + vertexCount]; % start counting face labels from the previous end label
        vertex = [vertex;jVertex]; % append new vertices
        if ~isequal(length(jColorData), size(jVertex,1))
            disp('here')
        end
        colorData = [colorData;jColorData];
    end
    % plot collision manifold
    %     patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
    patch('Faces', face, 'Vertices', vertex, 'FaceColor', [0,1,0], 'EdgeColor', 'none', 'FaceAlpha', .5);
    
end





if exist('atlasL4B')
    % plot the small collision manifold in F0 coordinates
    s = linspace(-1,1,100);
    tFinal = atlasL4B.MaxTau;
    t = linspace(0, tFinal, 200);
    [S,T] = meshgrid(s,t);
    evalData = [reshape(S,[],1), reshape(T,[],1)];
    
    face = [];
    vertex = [];
    colorData = [];
    for jChart = atlasL4B.Chart(1:end)
        bdVertices = [jChart.SpatialSpan(1), jChart.TimeSpan(1); jChart.SpatialSpan(1), jChart.TimeSpan(2);... % add 4 corners of domain
            jChart.SpatialSpan(2), jChart.TimeSpan(1); jChart.SpatialSpan(2), jChart.TimeSpan(2)];
        data = [jChart.intersectdomain(evalData);bdVertices];
        
        % map evaluations to F0 coordinates
        if isequal(jChart.RegType, 1)
            [X1,P1,Y1,Q1,~] = jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true);
            F0Data = CRTBP2reg([X1,P1,Y1,Q1],mu,-1);
            jColorData = 0.1*ones(size(data,1),1);
        elseif isequal(jChart.RegType, 2)
            [X2,P2,Y2,Q2,~] = jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true);
            F0Data = CRTBP2reg([X2,P2,Y2,Q2],mu,-2);
            jColorData = 0.5*ones(size(data,1),1);
        else
            [X0, P0, Y0, Q0,~,~] = jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true);
            F0Data = [X0, P0, Y0, Q0];
            jColorData = 0.9*ones(size(data,1),1);
        end
        
        X0 = F0Data(:,1);
        P0 = F0Data(:,2);
        Y0 = F0Data(:,3);
        Q0 = F0Data(:,4);
        jVertex = [mid(X0),mid(Y0),mid(P0)];
        jFace = delaunay(data(:,1),data(:,2)); % index triples for delaunay triangulation in space-time.
        vertexCount = size(vertex,1);
        face = [face;jFace + vertexCount]; % start counting face labels from the previous end label
        vertex = [vertex;jVertex]; % append new vertices
        if ~isequal(length(jColorData), size(jVertex,1))
            disp('here')
        end
        colorData = [colorData;jColorData];
    end
    % plot collision manifold
    %     patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
    patch('Faces', face, 'Vertices', vertex, 'FaceColor', [1,0,0], 'EdgeColor', 'none', 'FaceAlpha',.5);
end
axis([mu-1-.5, mu+.5, -1, 1, -3, 3])
dealfig()