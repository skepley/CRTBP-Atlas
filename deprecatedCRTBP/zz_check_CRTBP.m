%CHECK_REGCRTBP - testing for RegCRTBPChart class

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Mar-2019; Last revision: 18-Mar-2019

whichPC = 'mac';
switch whichPC
    case 'mac' % macbook
        computerPath = '/Users/sk2011/';
    case 'lenovo' % ubuntu
        computerPath = '/home/shane/';
end
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
clear all
clc

% ================================================== VERIFY NONRIGOROUS INTEGRATOR WORKS ==================================================
basis = 'Taylor';
initialTime = 0;
mu = .25;
parameter = [mu;1-mu];
isValid = false;

%% set up some initial data and advect manifold
thetaMid = 0;
thetaRadius = .0001;
N = 20;
truncation = [50,N];
expCoeff = exp(1i*thetaMid)*((1i*thetaRadius).^(0:N-1)./factorial(0:N-1));
xInitial = zeros(1,N);
pInitial = real(expCoeff);
yInitial = zeros(1,N);
qInitial = imag(expCoeff);
r1Initial = [1, zeros(1,N-1)];
r2Initial = [1, zeros(1,N-1)];
initialData = [xInitial;pInitial;yInitial;qInitial;r1Initial;r2Initial];


%% compute atlas for initial manifold
tDirection = 1;
tau = .1;
finalTau = .3;
bd = CRTBPChart(initialData, basis, initialTime, truncation, parameter, 'InitialScaling', tDirection*tau, 'Boundary', true);
atlas = Atlas(bd, tau, finalTau, @boundarycheck, @advectioncheck);
while ~isempty(atlas.LeafStack)
    atlas.growboundary()
end


%% ================================================== PLOT WITH RK45 ==================================================
close all
orbitFigure = figure;
orbitAxis = gca;
hold on

tValue = linspace(0,tDirection*tau,10)';
tRelative = linspace(0, 1 ,100)';
sRelative = linspace(-1, 1, 100)';
[S,T] = meshgrid(sRelative, tRelative);
evalData = [reshape(S, [],1), reshape(T, [], 1)];

% plot orbits
tf = linspace(0,finalTau,102);
s = 0;
ob = mid(atlas.orbit(s, tf));
plot3(ob(:,1), ob(:,3), ob(:,2), 'k.', 'LineWidth',3)

% rk45
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);
regVF = @(t,x)rk45vectorfield(t,x,mu);
initCondition = ob(1,:);
[~,sol] = ode45(regVF,tf,initCondition,odeOptions);
plot3(sol(:,1),sol(:,3),sol(:,2),'b')
view(-50,50)
dealfig()
