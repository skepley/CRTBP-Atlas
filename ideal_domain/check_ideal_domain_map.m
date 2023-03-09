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
%   email: s.kepley@vu.nl
%   Date: 31-Mar-2020; Last revision: 31-Mar-2020

close all
clear all
load ideal_domain_boundary_C3_equalmass
mu = 1/2;


%% ================================================== TEST F0 IDEAL DOMAIN FUNCTIONS ==================================================

% set up ideal domain functions
inP01 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P01Boundary{j}(1,:), P01Boundary{j}(2,:)), 1:length(P01Boundary)));  % u0 is not f1 ideal 
inP02 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P02Boundary{j}(1,:), P02Boundary{j}(2,:)), 1:length(P02Boundary)));  % u0 is not f2 ideal 
inP10 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P10Boundary{j}(1,:), P10Boundary{j}(2,:)), 1:length(P10Boundary)));  % u1 is not f1 ideal 
inP20 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P20Boundary{j}(1,:), P20Boundary{j}(2,:)), 1:length(P20Boundary)));  % u2 is not f2 ideal 
inP12 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P12Boundary{j}(1,:), P12Boundary{j}(2,:)), 1:length(P12Boundary)));  % u1 is not f2 ideal 
inP21 = @(u)any(arrayfun(@(j)inpolygon(u(1), u(2), P21Boundary{j}(1,:), P21Boundary{j}(2,:)), 1:length(P21Boundary)));  % u2 is not f1 ideal 


% ideal domain map
% F0 stay in F0 if: u in P01 and u in P02
% F0 swap to F1 if: x < 0 and (u not in both P01 and P02)
% F0 swap to F2 if: x > 0 and (u not in both P01 and P02)
idealdomainF0 = @(u)(u(1) <= 0 && (~inP01(u) || ~inP02(u))) + 2*(u(1) > 0 && (~inP01(u) || ~inP02(u)));

%  test the function
close all
scatterDensity = 50;
U1 = linspace(-0.8,0.8,scatterDensity);
U2 = linspace(-0.9,0.9,scatterDensity);
[U1,U2] = meshgrid(U1,U2);
UU = [reshape(U1, [], 1), reshape(U2, [], 1)];
% VV = arrayfun(@(j)idealdomainF0(UU(j,:)), 1:size(UU,1));
VV = arrayfun(@(j)equalmassidealdomain(UU(j,:), 0), 1:size(UU,1));

f0Ideal = [UU(VV==0,1), UU(VV==0,2)];  % set of f0 samples which should remmain in f0
f1Ideal = [UU(VV==1,1), UU(VV==1,2)];  % set of f0 samples which should map to f1
f2Ideal = [UU(VV==2,1), UU(VV==2,2)];  % set of f0 samples which should map to f1


% plot ideal point subsets in f0 coordinates
figure
hold on
scatter(f0Ideal(:,1), f0Ideal(:,2), 'r', 'filled') % f0 ideal
scatter(f1Ideal(:,1), f1Ideal(:,2), 'g', 'filled') % f1 ideal
scatter(f2Ideal(:,1), f2Ideal(:,2), 'b', 'filled') % f2 ideal

% % plot boundaries
% for j = 1:length(P01Boundary)
%     plot(P01Boundary{j}(1,:), P01Boundary{j}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)  % f0-f2 boundary top component
% end
% for j = 1:length(P02Boundary)
%     plot(P02Boundary{j}(1,:), P02Boundary{j}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)  % f0-f2 boundary top component
% end
title('Ideal domains in X0')


%% ================================================== TEST F1 IDEAL DOMAIN FUNCTIONS ==================================================
% map F0 ideal domain point subsets to f1 coordinates to get ground truth
w0 = sqrt(f0Ideal(:,1) + 1i*f0Ideal(:,2) - mu);
F1f0 = [real(w0), imag(w0)]; % image of f0Ideal under F1 maps to component of ideal domain in f1 coordinates which should be in f0 coordinates
w1 = sqrt(f1Ideal(:,1) + 1i*f1Ideal(:,2) - mu);
F1f1 = [real(w1), imag(w1)]; % image of f1Ideal under F1 maps to component of ideal domain in f1 coordinates which should remain in f1 coordinates
w2 = sqrt(f2Ideal(:,1) + 1i*f2Ideal(:,2) - mu);
F1f2 = [real(w2), imag(w2)]; % image of f0Ideal under F1 maps to component of ideal domain in f1 coordinates which should be in f0 coordinates
% 
% 
figure 
hold on
scatter(F1f0(:,1), F1f0(:,2), 'r', 'filled') % f0 ideal
scatter(F1f1(:,1), F1f1(:,2), 'g', 'filled') % f1 ideal
scatter(F1f2(:,1), F1f2(:,2), 'b', 'filled') % f2 ideal
% plot(P10Boundary{1}(1,:), P10Boundary{1}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
% plot(P10Boundary{2}(1,:), P10Boundary{2}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
% plot(P12Boundary{1}(1,:), P12Boundary{1}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)
% plot(P12Boundary{2}(1,:), P12Boundary{2}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)
% plot(y0f1(1,:), y0f1(2,:), 'k', 'LineWidth', 2.5)
% title('plot in X0 of ground truth')


% define boundary functions for the f1 ideal domains
% inF2strip = @(u)2*u(1).^2 > -mu + sqrt(mu^2 + 4*u(1).^2.*u(2).^2);  % returns true if point is inside the f2 ideal strip

% idealdomainF1 = @(u)(~inF2strip(u) && ~inP12(u)) + 2*(inF2strip(u) && ~inP10(u));

%  test the function
U1 = linspace(-0.8,0.8,scatterDensity);
U2 = linspace(-1.5,1.5,scatterDensity);
[U1,U2] = meshgrid(U1,U2);
UU = [reshape(U1, [], 1), reshape(U2, [], 1)];
% VV = arrayfun(@(j)idealdomainF1(UU(j,:)), 1:size(UU,1));
VV = arrayfun(@(j)equalmassidealdomain(UU(j,:), 1), 1:size(UU,1));


f0Ideal = [UU(VV==0,1), UU(VV==0,2)];  % set of f0 samples which should remmain in f0
f1Ideal = [UU(VV==1,1), UU(VV==1,2)];  % set of f0 samples which should map to f1
f2Ideal = [UU(VV==2,1), UU(VV==2,2)];  % set of f0 samples which should map to f1

% figure
% hold on
scatter(f0Ideal(:,1), f0Ideal(:,2), 'r', 'filled') % f0 ideal
scatter(f1Ideal(:,1), f1Ideal(:,2), 'g', 'filled') % f1 ideal
scatter(f2Ideal(:,1), f2Ideal(:,2), 'b', 'filled') % f2 ideal

% plot(P10Boundary{1}(1,:), P10Boundary{1}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
% plot(P10Boundary{2}(1,:), P10Boundary{2}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)
% plot(P12Boundary{1}(1,:), P12Boundary{1}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)
% plot(P12Boundary{2}(1,:), P12Boundary{2}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)
% plot(y0f1(1,:), y0f1(2,:), 'k', 'LineWidth', 2.5)
title('plot in X1 of my function')


%% set up F2 ideal domain map
% map F0 ideal domain point subsets to f2 coordinates to get ground truth

w0 = sqrt(f0Ideal(:,1) + 1i*f0Ideal(:,2) - mu + 1);
F2f0 = [real(w0), imag(w0); -real(w0), -imag(w0)]; % image of f0Ideal under F1 maps to component of ideal domain in f1 coordinates which should be in f0 coordinates
w1 = sqrt(f1Ideal(:,1) + 1i*f1Ideal(:,2) - mu + 1);
F2f1 = [real(w1), imag(w1); -real(w1), -imag(w1)]; % image of f1Ideal under F1 maps to component of ideal domain in f1 coordinates which should remain in f1 coordinates
w2 = sqrt(f2Ideal(:,1) + 1i*f2Ideal(:,2) - mu + 1);
F2f2 = [real(w2), imag(w2); -real(w2), -imag(w2)]; % image of f0Ideal under F1 maps to component of ideal domain in f1 coordinates which should be in f0 coordinates
% 
% 
figure 
hold on
scatter(F2f0(:,1), F2f0(:,2), 'r*') % f0 ideal
scatter(F2f1(:,1), F2f1(:,2), 'g*') % f1 ideal
scatter(F2f2(:,1), F2f2(:,2), 'b*') % f2 ideal
% 
% % plot boundaries
% for j = 1:length(P20Boundary)
%     plot(P20Boundary{j}(1,:), P20Boundary{j}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)  % f0-f2 boundary top component
% end
% for j = 1:length(P21Boundary)
%     plot(P21Boundary{j}(1,:), P21Boundary{j}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)  % f0-f2 boundary top component
% end


% 
% % set up f2 ideal domain maps
% inF1strip = @(u)2*u(1).^2 < 1 - mu + sqrt((1-mu)^2 + 4*u(1).^2.*u(2).^2);  % returns true if point is inside the f2 ideal strip
% idealdomainF2 = @(u)2*(~inF1strip(u) && ~inP21(u)) + (inF1strip(u) && ~inP20(u));


%  test the function
U1 = linspace(-1.5,1.5,scatterDensity);
U2 = linspace(-0.8,0.8,scatterDensity);
[U1,U2] = meshgrid(U1,U2);
UU = [reshape(U1, [], 1), reshape(U2, [], 1)];
% VV = arrayfun(@(j)idealdomainF2(UU(j,:)), 1:size(UU,1));
VV = arrayfun(@(j)equalmassidealdomain(UU(j,:), 2), 1:size(UU,1));



% f1f2B = @(x,y)x.^2 - 0.5*(1 - mu + sqrt((1-mu)^2 + 4*x.^2.*y.^2));
% Z = f1f2B(U1, U2);
% contour(U1, U2, Z, [0,0], 'r', 'LineWidth', 3)



f0Ideal = [UU(VV==0,1), UU(VV==0,2)];  % set of f0 samples which should remmain in f0
f1Ideal = [UU(VV==1,1), UU(VV==1,2)];  % set of f0 samples which should map to f1
f2Ideal = [UU(VV==2,1), UU(VV==2,2)];  % set of f0 samples which should map to f2

% figure
% hold on

% plot domains
scatter(f0Ideal(:,1), f0Ideal(:,2), 'r*') % f0 ideal
scatter(f1Ideal(:,1), f1Ideal(:,2), 'g*') % f1 ideal
scatter(f2Ideal(:,1), f2Ideal(:,2), 'b*') % f2 ideal
% scatter(-f0Ideal(:,1), -f0Ideal(:,2), 'r*') % f0 ideal mirrored
% 
% % plot boundaries
% for j = 1:length(P20Boundary)
%     plot(P20Boundary{j}(1,:), P20Boundary{j}(2,:), 'Color',[1,1,0], 'LineWidth', 2.5)  % f0-f2 boundary top component
% end
% for j = 1:length(P21Boundary)
%     plot(P21Boundary{j}(1,:), P21Boundary{j}(2,:), 'Color',[1,0,1], 'LineWidth', 2.5)  % f0-f2 boundary top component
% end
% plot(y0f2(1,:), y0f2(2,:), 'k', 'LineWidth', 2.5)

dealfig()
