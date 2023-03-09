%SAMPLE_IDEAL_DOMAIN - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       SAMPLE_IDEAL_DOMAIN description
%
%   Output:
%       SAMPLE_IDEAL_DOMAIN output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 14-Feb-2020; Last revision: 14-Feb-2020

%% ================================================== Define H0 energy gradient ==================================================
mu = 1/81;
r1 = @(u)sqrt((u(1) - mu).^2 + u(3).^2);
r2 = @(u)sqrt((u(1) + 1 - mu).^2 + u(3).^2);
gradH0 = @(u)2*[mu*(u(1) + 1 - mu).*(1 - r2(u)^(-3/2)) + (1 - mu).*(u(1) - mu).*(1 - r1(u).^(-3.2)),...
            -u(2),...
            mu*u(3).*(1 - r2(u)^(-3/2)) + (1 - mu).*u(3).*(1 - r1(u).^(-3.2)),...
            -u(4)];
% gradH0Norm = @(u)4*(u(2).^2 + u(4).^2 + mu^2*(1-2*r2(u).^(1/2) + r2(u).^(-1)) + (1-mu).^2*(1 - 2*r1(u).^(1/2) + r1(u).^(-1)));
gradH0Norm = @(u)norm(gradH0(u));


%% ================================================== Define D_u1(u0) ==================================================
N = @(u)u(3).^2 - u(1).^2;
D = @(u)u(1).^2 + u(3).^2;
Du0 = @(u)[2*u(1), 0, -2*u(3), 0;...
           u(2).*N(u)./(2*D(u).^2), u(1)./(2*D(u)), u(4).*N(u)./(2*D(u).^2), -u(3)./(2*D(u));...
           2*u(3), 0, 2*u(1), 0;...
           u(4).*N(u)./(2*D(u).^2), u(3)./(2*D(u)), -u(2).*N(u)./(2*D(u).^2), u(1)./(2*D(u))];


F1Full = @(u)CRTBP2reg(u, mu, 1);
projF1 = @(u6)u6(1:4); % project vectors in R^6 onto first 4 coordinates
F1 = @(u)projF1(F1Full(u)); 
gradH0RegNorm = @(u)norm(gradH0(u)*Du0(u));
% gradH0NormDiff = @(u)gradH0Norm(u) - gradH0RegNorm(u);
gradH0NormDiff = @(u)gradH0Norm(u)./gradH0RegNorm(u);



%% ================================================== Sample energy loss on configuration space ==================================================
U = @(xy)mu*(0.5*r2([xy(1), 0, xy(2), 0]).^2 + 1./r2([xy(1), 0, xy(2), 0])) + (1-mu)*(0.5*r1([xy(1), 0, xy(2), 0]).^2 + 1./r1([xy(1), 0, xy(2), 0])); % F0 potential
C = 3; % Fix an energy level to sample in
nSample = 100;
[X, Y] = meshgrid(linspace(-1,1,nSample), linspace(-1,1,nSample));
% Y = linspace(-1,1,nSample)';
velocityRadius = @(xy)sqrt(2*U(xy) - C);

% % SINGLE CHOICE OF VELOCITY AT THETA = PI/4
% P = @(xy)cos(pi/4)*velocityRadius(xy);
% Q = @(xy)sin(pi/4)*velocityRadius(xy);
% Z = arrayfun(@(x,y)gradH0NormDiff([x,P([x,y]),y,Q([x,y])]),X, Y); 

% AVERAGE OVER THE CIRCLE OF VELOCITIES IN THIS ENERGY
nVelocitySample = 100;
S1Sample = linspace(0,2*pi,nVelocitySample + 1);
S1Sample = S1Sample(1:end-1);
P = @(xy, theta)cos(theta)*velocityRadius(xy);
Q = @(xy, theta)sin(theta)*velocityRadius(xy);
gradH0MeanNormDiff = @(x,y)max(arrayfun(@(theta)gradH0NormDiff([x, P([x,y], theta), y, Q([x,y], theta)]), S1Sample));
Z = arrayfun(@(x,y)gradH0MeanNormDiff(x,y),X, Y); 

%%
close all
figure
hold on

minZ = min(Z(:));
maxZ = max(Z(:)); 
Zp = max(Z, -10);
Zp = min(Zp, 10);
contourf(X,Y,Zp,20)

colorbar
scatter([mu, mu-1],[0,0],'r*')
dealfig()
