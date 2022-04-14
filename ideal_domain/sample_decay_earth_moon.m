%SAMPLE_DECAY_EARTH_MOON - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       SAMPLE_DECAY_EARTH_MOON description
%
%   Output:
%       SAMPLE_DECAY_EARTH_MOON output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Apr-2020; Last revision: 18-Apr-2020

%% ================================================== SECTION 1 ==================================================

clear all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];
%% ================================================== SET PARAMETERS ==================================================

% Integrator and collision manifold parameters
mu = 1/81;  % earth-moon parameter value
C = 3.672580557517796329136672284221;  % L4 energy
M = 30; % time truncation
parameter = mu; 

% set up sampling functions
r1 = @(u)sqrt((u(1) - mu).^2 + u(3).^2);
r2 = @(u)sqrt((u(1) + 1 - mu).^2 + u(3).^2);
U = @(xy)mu*(0.5*r2([xy(1), 0, xy(2), 0]).^2 + 1./r2([xy(1), 0, xy(2), 0])) + (1-mu)*(0.5*r1([xy(1), 0, xy(2), 0]).^2 + 1./r1([xy(1), 0, xy(2), 0])); % F0 potential
sampleDomain = [-1.75, .9, 0, 1.62];
dx = 2/99;  % grid width to sample on
sampleX = sampleDomain(1):dx:sampleDomain(2);
sampleY = sampleDomain(3):dx:sampleDomain(4);
nSample = length(sampleX)*length(sampleY);
[X, Y] = meshgrid(sampleX, sampleY);
velocityRadius = @(xy)sqrt(2*U(xy) - C);

% velocity sampling functions for averaging over the circle in velocity coodinates for a fixed energy.
nVelocitySample = 5;  % number of samples to average over in velocity space
S1Sample = linspace(0,2*pi,nVelocitySample + 1);
S1Sample = S1Sample(1:end-1);
P = @(xy, theta)cos(theta)*velocityRadius(xy);
Q = @(xy, theta)sin(theta)*velocityRadius(xy);
%% ================================================== COMPUTE F0 DECAY MAP ==================================================
% load F0decay_C3_equalmass
% load F1decay_C3_equalmass

regType = 0;
meandecaymap = @(x,y)mean(arrayfun(@(theta)scalardecaymap([x; P([x,y], theta); y; Q([x,y], theta)], M, parameter, regType), S1Sample));
F0DecayData = arrayfun(@(x,y)meandecaymap(x,y), X, Y); 
save('F0decay_L4_earthmoon', 'F0DecayData');

%% ================================================== COMPUTE F1 DECAY MAP ==================================================
regType = 1;
parameter = [mu, C];
meandecaymap = @(x,y)mean(arrayfun(@(theta)scalardecaymap([x; P([x,y], theta); y; Q([x,y], theta)], M, parameter, regType), S1Sample));
F1DecayData = arrayfun(@(x,y)meandecaymap(x,y), X, Y); 
save('F1decay_L4_earthmoon', 'F1DecayData');

%% ================================================== COMPUTE F2 DECAY MAP ==================================================
regType = 2;
parameter = [mu, C];
meandecaymap = @(x,y)mean(arrayfun(@(theta)scalardecaymap([x; P([x,y], theta); y; Q([x,y], theta)], M, parameter, regType), S1Sample));
F2DecayData = arrayfun(@(x,y)meandecaymap(x,y), X, Y); 
save('F2decay_L4_earthmoon', 'F2DecayData');



%% ================================================== FILL IN DATA VIA REFLECTIONS ==================================================
Yr = [flipud(-Y(2:end,:));Y];
Xr = [X(2:end,:);X];
F1r = [flipud(F1DecayData(2:end,:)); F1DecayData];
F0r = [flipud(F0DecayData(2:end,:)); F0DecayData];
F2r = [flipud(F2DecayData(2:end,:)); F2DecayData];
% return
%% ================================================== GET VERTICES ALONG BREAKEVEN BOUNDARIES ==================================================

Z1 = F1r - F0r;
Z2 = F2r - F0r;
X = Xr;
Y = Yr;

% construct the boundary of P01 and P10 polygons
zeroLevelSet = contour(X,Y,Z1,[0,0],'r', 'LineWidth', 3); % use contour map to interpolate zero level set
close(gcf)
% initialize arrays of component polygons
P01Boundary = cell(1); % f1-f0 boundary in f0-coordinates
P10Boundary = cell(1); % f1-f0 boundary in f1-coordinates
P21Boundary = cell(1); % f1-f0 boundary in f2-coordinates

nComponent = 1;
while size(zeroLevelSet, 2) > 1
    componentIdx = zeroLevelSet(2,1) + 1;  % get final index of this component of zero level set
    P01Boundary{nComponent} = zeroLevelSet(:, 2:componentIdx);  % get vertices of this connected component of zero level set
    
    % apply F1 regularization map to find f1-f0 boundary in f1-coordinates
    x0 = P01Boundary{nComponent}(1,:);
    y0 = P01Boundary{nComponent}(2,:);
    w = sqrt(x0 + 1i*y0 - mu);
    P10Boundary{2*nComponent-1} = [real(w); imag(w)];  % positive branch of sqrt
    P10Boundary{2*nComponent} = [-real(w); -imag(w)];  % negative branch of sqrt

    
    % apply F2 regularization map to find f1-f0 boundary in f2-coordinates
    w = sqrt(x0 + 1i*y0 - (mu - 1));
    P21Boundary{2*nComponent-1} = [real(w); imag(w)];  % positive branch of sqrt
    P21Boundary{2*nComponent} = [-real(w); -imag(w)];  % positive branch of sqrt

    
    % increment to next component
    zeroLevelSet = zeroLevelSet(:, componentIdx + 1:end);  % remove this component from contour matrix
    nComponent = nComponent + 1;  % increment component
end

% construct the boundary of P02 and P20 polygons
zeroLevelSet = contour(X,Y,Z2,[0,0],'r', 'LineWidth', 3); % use contour map to interpolate zero level set
close(gcf)

% initialize arrays of component polygons
P02Boundary = cell(1); % f2-f0 boundary in f0-coordinates
P12Boundary = cell(1); % f2-f0 boundary in f1-coordinates
P20Boundary = cell(1); % f2-f0 boundary in f2-coordinates

nComponent = 1;
while size(zeroLevelSet, 2) > 1
    componentIdx = zeroLevelSet(2,1) + 1;  % get final index of this component of zero level set
    P02Boundary{nComponent} = zeroLevelSet(:, 2:componentIdx);  % get vertices of this connected component of zero level set
    
    % apply F2 regularization map to find f2-f0 boundary in f2-coordinates
    x0 = P02Boundary{nComponent}(1,:);
    y0 = P02Boundary{nComponent}(2,:);
    w = sqrt(x0 + 1i*y0 - (mu - 1));
    P20Boundary{2*nComponent-1} = [real(w); imag(w)];  % positive branch of sqrt
    P20Boundary{2*nComponent} = [-real(w); -imag(w)];  % negative branch of sqrt

    
    % apply F1 regularization map to find f2-f0 boundary in f1-coordinates
    w = sqrt(x0 + 1i*y0 - mu);
    P12Boundary{2*nComponent-1} = [real(w); imag(w)];  % positive branch of sqrt
    P12Boundary{2*nComponent} = [-real(w); -imag(w)];  % negtative branch of sqrt

    % increment to next component
    zeroLevelSet = zeroLevelSet(:, componentIdx + 1:end);  % remove this component from contour matrix
    nComponent = nComponent + 1;  % increment component
end

% construct the f1-f2 boundary which is the y-axis
y0axis = [zeros(1,100); linspace(-sampleDomain(4), sampleDomain(4), 100)];
w1 = sqrt(y0axis(1,:) + 1i*y0axis(2,:) - mu);
y0f1 = [real(w1); imag(w1)]; % image of y0 axis under F1
w2 = sqrt(y0axis(1,:) + 1i*y0axis(2,:) - (mu - 1));
y0f2 = [real(w2); imag(w2)]; % image of y0 axis under F2

return


% save('ideal_domain_boundary_C3_equalmass', 'P01Boundary', 'P10Boundary', 'P21Boundary', 'P02Boundary', 'P20Boundary', 'P12Boundary', 'y0f1', 'y0f2');




return
%% ================================================== PLOT DECAY DIFFERENCE ==================================================

% F10 CONTOUR LINES
Z1 = F1r - F0r;
X = Xr;
Y = Yr;
 
close all
figure
hold on

% plot contour lines
contourf(X,Y,Z1,25)

% plot zero level set
zeroLevelSet = contour(X,Y,Z1,[0,0],'r', 'LineWidth', 3);
component = cell(1);  % initialize array of component polygons
nComponent = 1;
while size(zeroLevelSet, 2) > 1
    componentIdx = zeroLevelSet(2,1) + 1;  % get final index of this component of zero level set
    component{nComponent} = zeroLevelSet(:, 2:componentIdx);
%     scatter(zeroPolyVert(1,:), zeroPolyVert(2,:), 'y')
    zeroLevelSet = zeroLevelSet(:, componentIdx + 1:end);  % remove this component from contour matrix
    nComponent = nComponent + 1;  % increment component
end

scatter([mu, mu-1], [0,0], 200, 'r', 'filled')  % plot primaries
colorbar
title('F1 - F0 decay')


% map ideal boundary to f1 coordinates 
figure
hold on
for j = 1:length(component)
    x0 = component{j}(1,:);
    y0 = component{j}(2,:);
    w = sqrt(x0 + 1i*y0 - mu);
    x1 = real(w);
    y1 = imag(w);
    plot(x1,y1);
end
dealfig()

return
%%  F20 CONTOUR LINES
Z2 = F2r - F0r;
figure
hold on

% plot contour lines
contourf(X,Y,Z2,25)

% plot zero level set
zeroLevelSet = contour(X,Y,Z2,[0,0],'r', 'LineWidth', 3);
zeroLevelSetComponents = cell(1);  % initialize array of component polygons
nComponent = 1;
while size(zeroLevelSet, 2) > 1
    componentIdx = zeroLevelSet(2,1) + 1;  % get final index of this component of zero level set
    component{nComponent} = zeroLevelSet(:, 2:componentIdx);
    P01Boundary = zeroLevelSet(:, 2:componentIdx);  % get vertices of this connected component of zero level set
%     scatter(zeroPolyVert(1,:), zeroPolyVert(2,:), 'y')
    zeroLevelSet = zeroLevelSet(:, componentIdx + 1:end);  % remove this component from contour matrix
    nComponent = nComponent + 1;  % increment component
end

scatter([mu, mu-1], [0,0], 200, 'r', 'filled')  % plot primaries
colorbar
title('F2 - F0 decay')

tCircle = linspace(0,2*pi, 100);
tCircle = tCircle(1:end-1);
xCircle = mu -1 + 0.5*cos(tCircle);
yCircle = 0.5*sin(tCircle);
plot(xCircle, yCircle, 'r--','LineWidth', 3)

% F12 CONTOUR LINES
Z12 = F1r - F2r;
figure
hold on

% plot contour lines
contourf(X,Y,Z12,25)

% plot zero level set
zeroLevelSet = contour(X,Y,Z12,[0,0],'r', 'LineWidth', 3);
zeroLevelSetComponents = cell(1);  % initialize array of component polygons
nComponent = 1;
while size(zeroLevelSet, 2) > 1
    componentIdx = zeroLevelSet(2,1) + 1;  % get final index of this component of zero level set
    component{nComponent} = zeroLevelSet(:, 2:componentIdx);
    zeroPolyVert = zeroLevelSet(:, 2:componentIdx);  % get vertices of this connected component of zero level set
%     scatter(zeroPolyVert(1,:), zeroPolyVert(2,:), 'y')
    zeroLevelSet = zeroLevelSet(:, componentIdx + 1:end);  % remove this component from contour matrix
    nComponent = nComponent + 1;  % increment component
end

scatter([mu, mu-1], [0,0], 200, 'r', 'filled')  % plot primaries
colorbar
title('F1 - F2 decay')
dealfig()

return
%% plot decay data
figure
contourf(X,Y,F0r,20)
colorbar
title('F0 decay')

figure
contourf(X,Y,F1r,20)
colorbar
title('F1 decay')

figure
contourf(X,Y,F2r,20)
colorbar
title('F2 decay')
dealfig()


return
%% define a function to check if a point is inside the ideal domain
isideal = @(x)inpolygon(x(1), x(2), component{1}(1,:), component{1}(2,:)) || inpolygon(x(1), x(2), component{2}(1,:), component{2}(2,:));
isideal([0,0])







%% ================================================== SECTION 2 ==================================================




