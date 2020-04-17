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
addpath(genpath([computerPath,'Dropbox/Matlab/SeqDESolver']))
addpath(genpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator/regCRTBP'))
addpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator')
rmpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator/CRTBP')
clear all
clc

%% ================================================== MAKE SOME CHOICES ==================================================
basis = 'Taylor';
initialTime = 0;
mu = .25; % small mass primary
isValid = false;
N = 10;
truncation = [25,N];
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);

% Integrator time stepping
tDirection = 1;
tau = .1;
finalTau = .3;
tf = linspace(0, finalTau, 250); % plotting time points


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

% Integrate multiple timesteps
bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tDirection*tau, 'boundary', true);
finalTau = .3;
atlas0 = Atlas(bd, tau, finalTau, @boundarycheck, @advectioncheck);
while ~isempty(atlas0.LeafStack)
    atlas0.growboundary()
end

% plot orbits
figure;
hold on
s = 0;
try
    ob = cell2mat(orbit.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']));
catch
    ob = mid(atlas0.orbit(s, tf));
end
plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',2)

% set up rk45 evaluation
VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
initCondition = ob(1,:)';
[~,sol] = ode45(VF,tf,initCondition,odeOptions);
plot3(sol(:,1),sol(:,3),sol(:,2),'go')
view(-50,50)


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
finalTau = .3;
atlas1 = Atlas(bd, tau, finalTau, @boundarycheck, @advectioncheck);
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
finalTau = .3;
atlas2 = Atlas(bd, tau, finalTau, @boundarycheck, @advectioncheck);
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

return
%% compute atlas for initial manifold
% change to interval arithmetic for testing rigorous solver
% d1 = sum(abs(ob(:,5) - 1./sqrt((ob(:,1) - mu).^2 + ob(:,3).^2)))
% d2 = sum(abs(ob(:,6) - 1./sqrt((ob(:,1)+1-mu).^2 + ob(:,3).^2)))

% % % %======== SINGLE ORBIT SEGMENT ===========
% % % finalTau = tau;
% % % orbit = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tDirection*tau);
% % % 
% % % bd = orbit.InitialData;
% % % x = bd(1);
% % % p = bd(2);
% % % y = bd(3);
% % % q = bd(4);
% % % R1 = bd(5);
% % % % R2 = bd(6);
% % % 
% % % x = orbit.Coordinate(1);
% % % p = orbit.Coordinate(2);
% % % y = orbit.Coordinate(3);
% % % q = orbit.Coordinate(4);
% % % R1 = orbit.Coordinate(5);
% % % % R2 = orbit.Coordinate(6);
% % % 
% % % 
% % % % u = (x+1-mu)*(x+1-mu) + y*y; % f0
% % % xx = x*x;
% % % yy = y*y;
% % %  u = (xx + yy)*(xx + yy) + 1 + 2*(xx - yy); % f1
% % % sq = sqrt(u);
% % % invsq = inv(sq);
% % % chk1 = sq*sq;
% % % sum(sum(abs(chk1.Coefficient - u.Coefficient)))
% % % chk2 = invsq*sq;
% % % log10(abs(chk2.Coefficient));
% % % sum(sum(abs(invsq.Coefficient-R1.Coefficient)))






%% ================================================== PLOT WITH RK45 ==================================================
% close all
% figure;
hold on

% plot orbits
tf = linspace(0, finalTau, 250);
s = 0;
try
    ob = cell2mat(orbit.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']));
catch
    ob = mid(atlas.orbit(s, tf));
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

return
%% compare with old integrator

addpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator/CRTBP')
tDirection = 1;
tau = .1;
obNew = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tDirection*tau);
obOld = CRTBPChart(initialData, basis, initialTime, truncation, [parameter,1-parameter], 'InitialScaling', tDirection*tau);

close all
figure
hold on
% 
% tValue = linspace(0,tDirection*tau,10)';
% tRelative = linspace(0, 1 ,10)';
% sRelative = linspace(-1, 1, 1)';
% [S,T] = meshgrid(sRelative, tRelative);
% evalData = [reshape(S, [],1), reshape(T, [], 1)];

% plot orbits
tf = linspace(0, tau, 102)';
s = 0;
evalData = [s*ones(size(tf)),tf];
[x0,p0,y0,~,~,~] = obOld.eval(evalData);
[x1,p1,y1,~,~,~] = obNew.eval(evalData);
plot3(x1, p1, y1, 'b', 'LineWidth',1)
plot3(x0, p0, y0, 'k.', 'LineWidth',3)
dealfig()

return
%% rk45
figure 
hold on
VFOld = @(t,x)rk45vectorfield(t,x,mu);
VFNew = @(t,x)rk45regvectorfield(t,x,mu,0);

% initCondition = ob(1,:)';
[~,solOld] = ode45(VFOld,tf,initCondition,odeOptions);
plot3(solOld(:,1),solOld(:,3),solOld(:,2),'r')
[~,solNew] = ode45(VFNew,tf,initCondition,odeOptions);
plot3(solNew(:,1),solNew(:,3),solNew(:,2),'b.')

solOld(1:5,:)
solNew(1:5,:)
sol(1:5,:)




return
