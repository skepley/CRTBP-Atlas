%CHECK_CONJUGACY - Verify the regularization formulae are correct
%
%   Description:
%       CHECK_CONJUGACY description
%
%   Output:
%       CHECK_CONJUGACY output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 04-Apr-2019; Last revision: 04-Apr-2019

addpath(genpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator'))
%% ================================================== CHECK BOTH RK45 FORMULAE ==================================================
close all
clear all
figure;
gca;
hold on
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);

%% MAKE SOME CHOICES
mu = 1/2; % small mass
mu1 = 1-mu;
regType = 1;

% F0 initial data
% x0i = 1;
% p0i = 2; 
% y0i = 3;
% q0i = 4;
x0i = rand(1,1);
p0i = rand(1,1);
y0i = rand(1,1);
q0i = rand(1,1);
r0i = 1./sqrt((x0i-mu).^2 + y0i.^2);
s0i = 1./sqrt((x0i+mu1).^2 + y0i.^2);
F0InitialData = [x0i;p0i;y0i;q0i;r0i;s0i];

%% Map to F1 then flow 
regInitCondition = CRTBP2reg(F0InitialData(1:4), mu, regType);
x10 = regInitCondition(1);
p10 = regInitCondition(2);
y10 = regInitCondition(3);
q10 = regInitCondition(4);
C = regInitCondition(6); % regularization energy
F1 = @(t,x)rk45regvectorfield(t,x,[mu,C],regType);

% integrate regularized vf 
F1Tau = 0.25;
[tau,F1sol] = ode45(F1, [0,F1Tau], regInitCondition(1:5)', odeOptions);

x1 = F1sol(:,1);
p1 = F1sol(:,2);
y1 = F1sol(:,3);
q1 = F1sol(:,4);

% plot3(x1(1:50:end), y1(1:50:end), p1(1:50:end), 'k', 'LineWidth', 2.5, 'LineStyle', ':')
F0Target = CRTBP2reg([x1(end), p1(end), y1(end), q1(end)], mu, -1);
scatter3(F0Target(1), F0Target(3), F0Target(2), 20, 'r*')
scatter3(x0i, y0i, p0i, 20, 'bo')


%% Flow in F0 then map to F1 
w = x1.^2 + y1.^2;
F0Tau = 4*trapz(tau,w);
tSpan = [0, F0Tau];

%  integrate F0
F0 = @(t,x)rk45regvectorfield(t,x,mu,0);
% [~,F0sol] = ode45(F0, tSpan, F0InitialData, odeOptions);
plot_orbit(F0, tSpan, F0InitialData, [1,3,2])


% % map orbit to F1
% x0 = F0sol(:,1);
% p0 = F0sol(:,2);
% y0 = F0sol(:,3);
% q0 = F0sol(:,4);
% F01Sol = CRTBP2reg([x0,p0,y0,q0],mu,1); % solutions integrated in F0 and mapped to F1
% x1 = F01Sol(:,1);
% y1 = F01Sol(:,3);
% p1 = F01Sol(:,2);
% plot3(x1, y1, p1, 'b', 'LineWidth', 1)

return

%% Map back to F0 and plot
figure
hold on
F0sol = CRTBP2reg([x1,p1,y1,q1],mu,-1);
xx0 = F0sol(:,1);
pp0 = F0sol(:,2);
yy0 = F0sol(:,3);
qq0 = F0sol(:,4);
plot3(xx0, yy0, pp0, 'b', 'LineWidth', 1)
plot3(x0i(1:50:end), y0i(1:50:end), p0(1:50:end), 'k', 'LineWidth', 2.5, 'LineStyle', ':')


dealfig()
axis tight
