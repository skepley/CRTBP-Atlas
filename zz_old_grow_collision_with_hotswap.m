%GROW_COLLISION_WITH_HOTSWAP - Compute a big collision manifold with dynamic hotswapping between F0, F1, F2.

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 21-May-2019; Last revision: 21-May-2019

clear all
close all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/SeqDESolver']))
addpath(genpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];

%% ================================================== MAKE SOME CHOICES ==================================================
tic
% Integrator and collision manifold parameters
basis = 'Taylor';
initialTime = 0; % always 0 for autonomous system
mu = 1/2; % small mass primary
isValid = false; % validate computations
N = 20; % spatial truncation
truncation = [40, N]; % time and space truncation spaces
C = 3.7;
parameter = [mu,C];


% Integrate the collision manifold for the large primary backwards in time.
regType = 1;
initialData = IDset(1-mu, N ,false); % set up some initial data on a F1 circle
maxGeneration = 20;
tau = -.1;
maxTau = -.001;
bd1 = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
% atlas1 = Atlas(bd1, tau, @boundarycheck, @advectioncheck, 'MaxGeneration', maxGeneration, 'MaxTau', maxTau);
atlas1 = Atlas(bd1, tau, @boundarycheck, @advectioncheck, 'MaxTau', maxTau);
while ~isempty(atlas1.LeafStack)
    fprintf('%0.4f \n', min([atlas1.LeafStack.TimeSpan]))
    atlas1.growboundary('RegTime', @regtime)
end

% % Integrate the collision manifold for the large primary forwards in time.
% regType = 1;
% initialData = IDset(mu, N ,false); % set up some initial data on a F1 circle
% maxGeneration = 20;
% tau = .1;
% maxTau = 5;
% bd2 = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
% atlas2 = Atlas(bd2, tau, @boundarycheck, @advectioncheck, 'MaxTau', maxTau);
% while ~isempty(atlas2.LeafStack)
%     atlas2.growboundary('RegTime', @regtime)
% end


% % Integrate the collision manifold for the small primary forwards in time.
% regType = 2;
% initialData = IDset(mu, N ,false); % set up some initial data on a F1 circle
% maxGeneration = 20;
% tau = .1;
% maxTau = .3;
% bd2 = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tau, 'boundary', true);
% atlas2 = Atlas(bd2, tau, @boundarycheck, @advectioncheck, 'MaxTau', maxTau);
% while ~isempty(atlas2.LeafStack)
%     fprintf('%0.4f \n', max([atlas2.LeafStack.TimeSpan]))
%     atlas2.growboundary('RegTime', @regtime)
% end

toc
return
% save([savePath, 'collisions_3_3dot7_equalmass'])
%% PLOT THE MANIFOLDS
close all
figure
hold on

% plot primaries
primarySize = 500;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','b')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')

% plot bounding circles
bdtheta = linspace(0,2*pi,500);
xb = 0.25*cos(bdtheta);
yb = 0.25*sin(bdtheta);
pb = zeros(size(bdtheta));
xb1 = xb.^2 - yb.^2 + mu;
yb1 = 2*xb.*yb;
plot3(xb1,yb1,zeros(size(xb)),'y','LineWidth',3)
xb2 = xb.^2 - yb.^2 + (mu-1);
yb2 = 2*xb.*yb;
plot3(xb2,yb2,zeros(size(xb)),'y','LineWidth',3)


% plot the large collision manifold in F0 coordinates
s = linspace(-1,1,100);
% tFinal = max([atlas1.Chart.TimeSpan]);
tFinal = atlas1.MaxTau;
t = linspace(0, tFinal, 200);
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];


if exist('atlas1')
    
    face = [];
    vertex = [];
    colorData = [];
    
    for jChart = atlas1.Chart(1:end)
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
    patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
    % patch('Faces', face, 'Vertices', vertex, 'FaceColor', [1,0,0], 'EdgeColor', 'none', 'FaceAlpha',.5);
end


if exist('atlas2')
    % plot the small collision manifold in F0 coordinates
    s = linspace(-1,1,100);
    tFinal = atlas2.MaxTau;
    t = linspace(0, tFinal, 200);
    [S,T] = meshgrid(s,t);
    evalData = [reshape(S,[],1), reshape(T,[],1)];
    
    face = [];
    vertex = [];
    colorData = [];
    for jChart = atlas2.Chart(1:end)
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
    patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
    % patch('Faces', face, 'Vertices', vertex, 'FaceColor', [0,0,1], 'EdgeColor', 'none', 'FaceAlpha',.5);
end


axis([mu-1-.5, mu+.5, -1, 1, -10, 10])
dealfig()
