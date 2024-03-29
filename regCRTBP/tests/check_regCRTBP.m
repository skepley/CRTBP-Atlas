%CHECK_REGCRTBP - testing for RegCRTBPChart class

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 18-Mar-2019; Last revision: 22-Apr-2020

clear all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))

%% ================================================== MAKE SOME CHOICES ==================================================


% Integrator parameters
N = 10; % spatial truncation
M = 25; % temporal truncation
truncation = [M, N]; % set truncation vector
subDivideParameter = [4, .1, 4]; % default is [4, .1, 4]
basis = 'Taylor';
mu = .25; % small mass primary
isValid = false;
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);

% Integrator time stepping
initialTime = 0;
tDirection = 1;
tauGuess = 0.1; % initial timestep to scale each Taylor step
maxTau = .3;
tf = linspace(0, maxTau, 250); % plotting time points


%% ================================================== CHECK f_0 AGAINST RUNGE KUTTA ==================================================
regType = 0;
parameter = mu;
% set up some initial data on a line
p1 = [0.7104, 0.2973, 0.2554, 0.2068];
p2 = p1 + 1e-4*ones(size(p1));
initialData = cat(2, 0.5*[p1+p2; p2-p1]', zeros(4,N-2));
if isValid
    initialData = intval(initialData);
end
bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);

% Integrate multiple timesteps
bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau, 'SubDivideParameter', subDivideParameter);
A = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
while ~isempty(A.LeafStack)
    A.growboundary()
end

% plot orbits
figure;
hold on
s = 0;
try
    ob = cell2mat(orbit.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']));
catch
    ob = mid(A.orbit(s, tf));
end
plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',2)

% set up rk45 evaluation
VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
initCondition = ob(1,:)';
[~,sol] = ode45(VF,tf,initCondition,odeOptions);
plot3(sol(:,1),sol(:,3),sol(:,2),'go')
view(-50,50)

return
%% ================================================== CHECK f_1 AGAINST RUNGE KUTTA ==================================================
regType = 1;
C = 3;
parameter = [mu,C];

% set up some initial data on a circle
thetaMid = 0;
thetaRadius = .0001;
expCoeff = sqrt(8*mu)*exp(1i*thetaMid)*((1i*thetaRadius).^(0:N-1)./factorial(0:N-1));
xInitial = zeros(1,N);
pInitial = real(expCoeff);
yInitial = zeros(1,N);
qInitial = imag(expCoeff);
initialData = [xInitial;pInitial;yInitial;qInitial];
if isValid
    initialData = intval(initialData);
end

% Integrate multiple timesteps
bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tDirection*tau, 'boundary', true);
maxTau = .3;
atlas1 = Atlas(bd, tau, maxTau, @boundarycheck, @advectioncheck);
while ~isempty(atlas1.LeafStack)
    atlas1.growboundary()
end

% plot orbits
figure;
hold on
s = 0;
try
    ob = cell2mat(orbit.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']));
catch
    ob = mid(atlas1.orbit(s, tf));
end
plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',2)

% set up rk45 evaluation
VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
initCondition = ob(1,:)';
[~,sol] = ode45(VF,tf,initCondition,odeOptions);
plot3(sol(:,1),sol(:,3),sol(:,2),'go')
view(-50,50)

%% ================================================== CHECK f_2 AGAINST RUNGE KUTTA ==================================================
regType = 2;
C = 3;
parameter = [mu,C];

% set up some initial data on a circle
thetaMid = 0;
thetaRadius = .0001;
expCoeff = sqrt(8*mu)*exp(1i*thetaMid)*((1i*thetaRadius).^(0:N-1)./factorial(0:N-1));
xInitial = zeros(1,N);
pInitial = real(expCoeff);
yInitial = zeros(1,N);
qInitial = imag(expCoeff);
initialData = [xInitial;pInitial;yInitial;qInitial];
if isValid
    initialData = intval(initialData);
end

% Integrate multiple timesteps
bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tDirection*tau, 'boundary', true);
maxTau = .3;
atlas2 = Atlas(bd, tau, maxTau, @boundarycheck, @advectioncheck);
while ~isempty(atlas2.LeafStack)
    atlas2.growboundary()
end

% plot orbits
% figure;
hold on
s = 0;
try
    ob = cell2mat(orbit.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']));
catch
    ob = mid(atlas2.orbit(s, tf));
end
plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',2)

% set up rk45 evaluation
VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
initCondition = ob(1,:)';
[~,sol] = ode45(VF,tf,initCondition,odeOptions);
plot3(sol(:,1),sol(:,3),sol(:,2),'go')
view(-50,50)
dealfig()


%% ================================================== PLOT WITH RK45 ==================================================
% close all
% figure;
hold on

% plot orbits
tf = linspace(0, maxTau, 250);
s = 0;
try
    ob = cell2mat(orbit.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']));
catch
    ob = mid(A.orbit(s, tf));
end
plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',1)

% d1 = sum(abs(ob(:,5) - 1./sqrt((ob(:,1) - mu).^2 + ob(:,3).^2)))
% d2 = sum(abs(ob(:,6) - 1./sqrt((ob(:,1)+1-mu).^2 + ob(:,3).^2)))

% set up rk45 evaluation
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);
VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
initCondition = ob(1,:)';
[~,sol] = ode45(VF,tf,initCondition,odeOptions);
plot3(sol(:,1),sol(:,3),sol(:,2),'go')
view(-50,50)
dealfig()


return
% integrate with Jay's code
addpath('/Users/sk2011/Dropbox/TaylorDocs/CRFBP_regularization/codes')
% 
% tspan = linspace(0, 0.3, 5000);
% IC = [x,pz,y,qz]

%Integraton:
options=odeset('RelTol',1e-13,'AbsTol',1e-13);
[t,thisOrbit_reg] = ode45('crtb_Field', ...
    tf, initCondition(1:4), options, flag, 1-mu, mu);

plot3(thisOrbit_reg(:,1), thisOrbit_reg(:,3), thisOrbit_reg(:,2), 'go', 'LineWidth',1)

% scatter3([1-mu,mu],[0,0], [0,0], 100,'b','filled')

% crtb_Field(1,initCondition(1:4),options, flag, 1-mu,mu)

