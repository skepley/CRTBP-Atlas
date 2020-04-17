%CHECK_IDEAL_DOMAIN_MAP - Verify and test ideal domain code is correct
%
%   Description:
%       CHECK_IDEAL_DOMAIN_MAP description
%
%   Output:
%       CHECK_IDEAL_DOMAIN_MAP output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 31-Mar-2020; Last revision: 31-Mar-2020

close all
clear all
%% ================================================== TEST IDEAL DOMAIN FUNCTIONS ==================================================
load ideal_domain_boundary_C3_equalmass
mu = 1/2;

% set up ideal domain functions
inP01 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P01Boundary{j}(1,:), P01Boundary{j}(2,:)), 1:length(P01Boundary)));  % returns true if u in P01
inP02 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P02Boundary{j}(1,:), P02Boundary{j}(2,:)), 1:length(P02Boundary)));  % returns true if u in P02
inP10 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P10Boundary{j}(1,:), P10Boundary{j}(2,:)), 1:length(P10Boundary)));  % returns true if u in P10
inP20 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P20Boundary{j}(1,:), P20Boundary{j}(2,:)), 1:length(P20Boundary)));  % returns true if u in P20

% ================================================== ideal domain maps: ==================================================
% F0 stay in F0: u in P01 and u in P02
% F0 swap to F1: x < 0 and (u not in both P01 and P02)
% F0 swap to F2: x > 0 and (u not in both P01 and P02)
idealdomainF0 = @(u)(u(1) <= 0 && (~inP01(u) || ~inP02(u))) + 2*(u(1) > 0 && (~inP01(u) || ~inP02(u)));

% F1 stay in F1: u not in P10
% F1 swap to F0: u in P10
idealdomainF1 = @(u)~inP10(u);

% F2 stay in F2: u not in P20
% F2 swap to F0: u in P20
idealdomainF2 = @(u)2*~inP20(u);


%%  test
close all
U1 = linspace(-1,1,40);
U2 = linspace(-1.3,1.3,40);
[U1,U2] = meshgrid(U1,U2);
UU = [reshape(U1, [], 1), reshape(U2, [], 1)];
VV = arrayfun(@(j)idealdomainF0(UU(j,:)), 1:size(UU,1));

f0Ideal = [UU(VV==0,1), UU(VV==0,2)];  % set of f0 samples which should remmain in f0
f1Ideal = [UU(VV==1,1), UU(VV==1,2)];  % set of f0 samples which should map to f1
f2Ideal = [UU(VV==2,1), UU(VV==2,2)];  % set of f0 samples which should map to f1

% plot ideal point subsets in f0 coordinates
figure
hold on
plot(P01Boundary{1}(1,:), P01Boundary{1}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
plot(P01Boundary{2}(1,:), P01Boundary{2}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
plot(P02Boundary{1}(1,:), P02Boundary{1}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)
plot(P02Boundary{2}(1,:), P02Boundary{2}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)

scatter(f0Ideal(:,1), f0Ideal(:,2), 'r*') % f0 ideal
scatter(f1Ideal(:,1), f1Ideal(:,2), 'g*') % f1 ideal
scatter(f2Ideal(:,1), f2Ideal(:,2), 'b*') % f2 ideal


% map ideal point subsets to f1 coordinates
w0 = sqrt(f0Ideal(:,1) + 1i*f0Ideal(:,2) - mu);
F1f0 = [real(w0), imag(w0)]; % image of f0Ideal under F1 maps to component of ideal domain in f1 coordinates which should be in f0 coordinates
w1 = sqrt(f1Ideal(:,1) + 1i*f1Ideal(:,2) - mu);
F1f1 = [real(w1), imag(w1)]; % image of f1Ideal under F1 maps to component of ideal domain in f1 coordinates which should remain in f1 coordinates
w2 = sqrt(f2Ideal(:,1) + 1i*f2Ideal(:,2) - mu);
F1f2 = [real(w2), imag(w2)]; % image of f0Ideal under F1 maps to component of ideal domain in f1 coordinates which should be in f0 coordinates


figure
hold on
scatter(F1f0(:,1), F1f0(:,2), 'r*') % f0 ideal
scatter(F1f1(:,1), F1f1(:,2), 'g*') % f1 ideal
scatter(F1f2(:,1), F1f2(:,2), 'b*') % f2 ideal
plot(P10Boundary{1}(1,:), P10Boundary{1}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
plot(P10Boundary{2}(1,:), P10Boundary{2}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
plot(P12Boundary{1}(1,:), P12Boundary{1}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)
plot(P12Boundary{2}(1,:), P12Boundary{2}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)
% plot(y0f1(1,:), y0f1(2,:), 'k', 'LineWidth', 2.5)
scatter(y0f1(1,:), y0f1(2,:), 20, 'k', 'filled')



%% set up ideal domain maps
x = @(x1,y1)x1 - sqrt(0.5*(-mu + sqrt(mu^2 + 4*x1.^2.*y1.^2)));
X1 = linspace(0,1,100);
Y1 = linspace(-1.5,1.5,100);
[X1,Y1] = meshgrid(X1,Y1);
z = x(X1,Y1);
figure
hold on
contourf(X1,Y1,z,25)
contour(X1,Y1,z,[0,0],'r', 'LineWidth',3)
colorbar()

dealfig()

