%CHECK_HOTSWAP - verify the map between regularized and standard coordinates works for Chart objects. 
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 06-Apr-2019; Last revision: 06-Apr-2019

% ================================================== SECTION 1 ==================================================
clear all
whichPC = 'mac';
switch whichPC
    case 'mac' % macbook
        computerPath = '/Users/sk2011/';
    case 'lenovo' % ubuntu
        computerPath = '/home/shane/';
end
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator'))
clc
% close all
%% ================================================== MAKE SOME CHOICES ==================================================
basis = 'Taylor';
initialTime = 0;
mu = .25; % small mass primary
isValid = false;
N = 10;
truncation = [25,N];
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);

% Integrator time stepping
tau = .1;
finalTau = .4;
tf = linspace(0, finalTau, 250); % plotting time points

%% ================================================== SWAP FROM F_2 ---> F_0 ==================================================
swapFrom = 2;
swapTo = 0;
C = 3;
parameter = [mu,C];
initialData = IDset(mu, N ,false); % set up some initial data on a circle

% Integrate multiple timesteps
bd1 = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, swapFrom, 'InitialScaling', tau, 'boundary', true);
atlas = Atlas(bd1, tau, finalTau, @boundarycheck, @advectioncheck);
while ~isempty(atlas.LeafStack)
    atlas.growboundary()
end


%% ================================================== PLOT IN F_2 SPACE ==================================================
figure;
hold on
tValue = linspace(0,tau,50)';
tRelative = linspace(0, 1, 25)';
sRelative = linspace(-1, 1, 25)';
[S,T] = meshgrid(sRelative, tRelative);
evalData = [reshape(S, [],1), reshape(T, [], 1)];

% plot forward 
% interior and initial boundary
for j = 1:atlas.Size
    jChart = atlas.Chart(j);
    [x,P1,y,Q1,r] = jChart.eval(evalData); % Taylor integrator
    if jChart.Generation == 0
        [x0,p0,y0,q0,z0] = jChart.eval([sRelative, zeros(size(sRelative))]);
        plot3(mid(x0),mid(y0),mid(p0), 'k', 'LineWidth',3)
    end
    scatter3(mid(x),mid(y), mid(P1), 5,'filled', 'b')
end

% plot terminal boundaries
termChart = RegCRTBPChart;
for jChart = 1:length(atlas.CrashStack)
    parentChart = atlas.CrashStack(jChart).ParentHandle;
    termChart(jChart) = parentChart;
    [x1,p1,y1,q1,z1] = parentChart.eval([sRelative, ones(size(sRelative))]);
    plot3(mid(x1),mid(y1),mid(p1), 'r', 'LineWidth',3)
end

% plot bounding circle
bdtheta = linspace(0,2*pi,500);
xb = 0.5*cos(bdtheta);
yb = 0.5*sin(bdtheta);
pb = zeros(size(bdtheta));
plot3(xb,yb,pb,'y','LineWidth',3)

%%  ============================= PLOT FORWARD TIME F_2 COLLISION MANIFOLD IN F_0 COORDINATES ===============================
s = linspace(-1,1,100);
t = linspace(.1,finalTau,100);
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];

face = [];
vertex = [];
colorData = [];

for jChart = atlas.Chart
    bdVertices = [jChart.SpatialSpan(1), jChart.TimeSpan(1); jChart.SpatialSpan(1), jChart.TimeSpan(2);... % add 4 corners of domain
        jChart.SpatialSpan(2), jChart.TimeSpan(1); jChart.SpatialSpan(2), jChart.TimeSpan(2)];
    data = [jChart.intersectdomain(evalData);bdVertices];
    [X1,P1,Y1,Q1,~] = jChart.eval(data,'GlobalTime', true, 'GlobalSpace', true);
    F0Data = CRTBP2reg([X1,P1,Y1,Q1],mu,-2);
    X0 = F0Data(:,1);
    Y0 = F0Data(:,3);
    P0 = F0Data(:,2);
    jColorData = data(:,2);
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



% plot it
% close all
% figure
hold on

% plot primaries
primarySize = 500;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','b')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')

% plot collision manifold
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);

view(10,9)
dealfig() 
colorbar()
set(gca,'ZLim',[-15,15])


% plot bounding circle
bdtheta = linspace(0,2*pi,500);
xb = 0.5*cos(bdtheta);
yb = 0.5*sin(bdtheta);
pb = zeros(size(bdtheta));
xb0 = xb.^2 - yb.^2 + mu-1;
yb0 = 2*xb.*yb;
plot3(xb0,yb0,zeros(size(xb)),'y','LineWidth',3)


% map to new coordinates and plot
orbit = atlas.Chart(end);
bd = fixtime(orbit,1);
bd0 = bd.hotswap(0);



% plot initial data in F0 coordinates
[x,p,y,~,~,~] = bd0.eval(s'); % Taylor integrator
aa = atlas.Chart(end);
[xx,pp,yy,qq,~] = aa.eval([s',ones(length(s),1)]);
newData = CRTBP2reg([xx,pp,yy,qq],mu,-2);
X0 = newData(:,1);
Y0 = newData(:,3);
P1 = newData(:,2);
scatter3(mid(x),mid(y), mid(p), 5,'filled', 'b')
scatter3(X0,Y0,P1,5,'r')

%%
newFinalTau = .2;
atlas0 = Atlas(bd0, tau, newFinalTau, @boundarycheck, @advectioncheck);
while ~isempty(atlas0.LeafStack)
    atlas0.growboundary()
end

% plot patches for F0 atlas
s = linspace(-1,1,100);
t = linspace(0,newFinalTau,100);
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];

face = [];
vertex = [];
colorData = [];

for jChart = atlas0.Chart
    bdVertices = [jChart.SpatialSpan(1), jChart.TimeSpan(1); jChart.SpatialSpan(1), jChart.TimeSpan(2);... % add 4 corners of domain
        jChart.SpatialSpan(2), jChart.TimeSpan(1); jChart.SpatialSpan(2), jChart.TimeSpan(2)];
    data = [jChart.intersectdomain(evalData);bdVertices];
    [x,p,y,~,~,~] = jChart.eval(data,'GlobalTime',true,'GlobalSpace',true);
    jColorData = data(:,2);
    jVertex = [mid(x),mid(y),mid(p)];
    jFace = delaunay(data(:,1),data(:,2)); % index triples for delaunay triangulation in space-time.
    vertexCount = size(vertex,1);
    face = [face;jFace + vertexCount]; % start counting face labels from the previous end label
    vertex = [vertex;jVertex]; % append new vertices
    if ~isequal(length(jColorData), size(jVertex,1))
        disp('here')
    end
    colorData = [colorData;jColorData];
end

patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
return
%%  ============================= PLOT BACKWARD TIME F_2 COLLISION MANIFOLD FROM F_0 TO F_2 COORDINATES ===============================

% map to new coordinates and plot
orbit1 = atlas0.Chart(end);
bd = fixtime(orbit1,1);
bd2 = bd.hotswap(2);

% plot initial data in F1 coordinates
[x2,p2,y2,~,~,~] = bd2.eval(s'); 
[xx,pp,yy,qq,~,~] = orbit1.eval([s',ones(length(s),1)]);
newData = CRTBP2reg([xx,pp,yy,qq],mu,2);
X2 = newData(:,1);
P2 = newData(:,2);
Y2 = newData(:,3);
scatter3(mid(x2), mid(y2), mid(p2), 5,'filled', 'b')
scatter3(X2,Y2,P2,5,'r')

%  ============================= CONTINUE INTEGRATION OF COLLISION MANIFOLD IN F_1 COORDINATES ===============================
tRelative = linspace(0, -1, 25)';
sRelative = linspace(-1, 1, 25)';
[S,T] = meshgrid(sRelative, tRelative);
evalData = [reshape(S, [],1), reshape(T, [], 1)];

newFinalTau = -.3;
tau = -.1;
atlas2 = Atlas(bd2, tau, newFinalTau, @boundarycheck, @advectioncheck);
while ~isempty(atlas2.LeafStack)
    atlas2.growboundary()
end

% plot backward time manifolds in F2 
s = linspace(-1,1,100);
t = linspace(0,newFinalTau,100);
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];

face = [];
vertex = [];
colorData = [];

for jChart = atlas2.Chart
    bdVertices = [jChart.SpatialSpan(1), jChart.TimeSpan(1); jChart.SpatialSpan(1), jChart.TimeSpan(2);... % add 4 corners of domain
        jChart.SpatialSpan(2), jChart.TimeSpan(1); jChart.SpatialSpan(2), jChart.TimeSpan(2)];
    data = [jChart.intersectdomain(evalData);bdVertices];
    [x,p,y,~,~,~] = jChart.eval(data,'GlobalTime',true,'GlobalSpace',true);
    jColorData = data(:,2);
    jVertex = [mid(x),mid(y),mid(p)];
    jFace = delaunay(data(:,1),data(:,2)); % index triples for delaunay triangulation in space-time.
    vertexCount = size(vertex,1);
    face = [face;jFace + vertexCount]; % start counting face labels from the previous end label
    vertex = [vertex;jVertex]; % append new vertices
    if ~isequal(length(jColorData), size(jVertex,1))
        disp('here')
    end
    colorData = [colorData;jColorData];
end
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
dealfig()

