%SAMPLE_DECAY_MAP - Sample and evaluate the decay map for the 3 vector fields to determine the boundary of the ideal domain.
%
%   Description:
%       SAMPLE_DECAY_MAP Samples a map of the form (x0, y0) ---> tau which gives the natural timestep of the Taylor
%           integration of an IVP for f0 with initial condition (x0, y0, p0, q0). For a fixed energy level set, p0, q0
%           must be chosen to lie on a circle of fixed radius. This map is computed uniformly on this circle and averaged.
%
%   Output:
%       F0DecayData.mat, F1DecayData.mat
%
%   Other m-files required: @RegCRTBPChart
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 19-Feb-2020; Last revision: 9-Mar-2020


clear all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath('/Users/shane/Dropbox/Regularisation3bp/CRTBP Atlas'))
savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];
%% ================================================== SET PARAMETERS ==================================================

% Integrator and collision manifold parameters
mu = 1/2;
M = 30; % time truncation
C = 3; % Fix an energy level to sample in
parameter = mu;

% set up sampling functions
r1 = @(u)sqrt((u(1) - mu).^2 + u(3).^2);
r2 = @(u)sqrt((u(1) + 1 - mu).^2 + u(3).^2);
U = @(xy)mu*(0.5*r2([xy(1), 0, xy(2), 0]).^2 + 1./r2([xy(1), 0, xy(2), 0])) + (1-mu)*(0.5*r1([xy(1), 0, xy(2), 0]).^2 + 1./r1([xy(1), 0, xy(2), 0])); % F0 potential
sampleDomain = [-.9, .6, 0, 1.3];
dx = 2/99;  % grid width to sample on
sampleX = sampleDomain(1):dx:sampleDomain(2);
sampleY = sampleDomain(3):dx:sampleDomain(4);
nSample = length(sampleX)*length(sampleY);
[X, Y] = meshgrid(sampleX, sampleY);
velocityRadius = @(xy)sqrt(2*U(xy) - C);

% velocity sampling functions for averaging over the circle in velocity coodinates for a fixed energy.
nVelocitySample = 10;  % number of samples to average over in velocity space
S1Sample = linspace(0,2*pi,nVelocitySample + 1);
S1Sample = S1Sample(1:end-1);
P = @(xy, theta)cos(theta)*velocityRadius(xy);
Q = @(xy, theta)sin(theta)*velocityRadius(xy);

%% ================================================== COMPUTE F0 DECAY MAP ==================================================
tic 
regType = 0;
meandecaymap = @(x,y)mean(arrayfun(@(theta)scalardecaymap([x; P([x,y], theta); y; Q([x,y], theta)], M, parameter, regType), S1Sample));
F0DecayData = arrayfun(@(x,y)meandecaymap(x,y), X, Y); 
save('F0decay_C3_equalmass', 'F0DecayData');
toc

%% ================================================== COMPUTE F1 DECAY MAP ==================================================
regType = 1;
parameter = [mu, C];
meandecaymap = @(x,y)mean(arrayfun(@(theta)scalardecaymap([x; P([x,y], theta); y; Q([x,y], theta)], M, parameter, regType), S1Sample));
F1DecayData = arrayfun(@(x,y)meandecaymap(x,y), X, Y); 
save('F1decay_C3_equalmass', 'F1DecayData');

%% ================================================== COMPUTE F2 DECAY MAP ==================================================
regType = 2;
parameter = [mu, C];
meandecaymap = @(x,y)mean(arrayfun(@(theta)scalardecaymap([x; P([x,y], theta); y; Q([x,y], theta)], M, parameter, regType), S1Sample));
F1DecayData = arrayfun(@(x,y)meandecaymap(x,y), X, Y); 
save('F1decay_C3_equalmass', 'F1DecayData');

return
%% ================================================== PLOT DECAY DIFFERENCE ==================================================

% load F0DecayData
% load F1DecayData
Z = F1DecayData - F0DecayData;

close all
figure
hold on

% plot contour lines
contourf(X,Y,Z,25)

% plot zero level set
zeroLevelSet = contour(X,Y,Z,[0,0],'r', 'LineWidth', 3);
zeroLevelSetComponents = cell(1);  % initialize array of component polygons
nComponent = 1;
while size(zeroLevelSet, 2) > 1
    componentIdx = zeroLevelSet(2,1) + 1;  % get final index of this component of zero level set
    component{nComponent} = zeroLevelSet(:, 2:componentIdx);
    zeroPolyVert = zeroLevelSet(:, 2:componentIdx);  % get vertices of this connected component of zero level set
    scatter(zeroPolyVert(1,:), zeroPolyVert(2,:), 'y')
    zeroLevelSet = zeroLevelSet(:, componentIdx + 1:end);  % remove this component from contour matrix
    nComponent = nComponent + 1;  % increment component
end

scatter([mu, mu-1], [0,0], 200, 'r', 'filled')  % plot primaries
colorbar
title('F1 - F0 decay')

tCircle = linspace(0,2*pi, 100);
tCircle = tCircle(1:end-1);
xCircle = mu + 0.5*cos(tCircle);
yCircle = 0.5*sin(tCircle);
plot(xCircle, yCircle, 'r--','LineWidth', 3)

% define a function to check if a point is inside the ideal domain
isideal = @(x)inpolygon(x(1), x(2), component{1}(1,:), component{1}(2,:)) || inpolygon(x(1), x(2), component{2}(1,:), component{2}(2,:));
isideal([0,0])

% check mapping ideal boundary to f1 coordinates
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
% figure
% contourf(X,Y,F0DecayData,20)
% colorbar
% title('F0 decay')
% 
% figure
% contourf(X,Y,scalarDecayData,20)
% colorbar
% title('F1 decay')

% define contour polygon

dealfig()

