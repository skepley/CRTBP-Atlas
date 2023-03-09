%CHECK_SHOOTING_TEMPLATE -A bit of a playground for developing the shooting template class
%
%   Description:
%       CHECK_SHOOTING_TEMPLATE description
%
%   Output:
%       CHECK_SHOOTING_TEMPLATE output
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 03-Mar-2023;


%% Get some data to test with
clearvars
close all
computerPath = pcpath('mac');
addpath('/users/shane/dropbox/matlab/tools/exportfig')
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
load newl4_equalmass_local_manifolds % get local manifold data
if ~exist('connections', 'var')
    load chk_connection_save.mat
end
mu = connections(1).Parameter;
C = connections(45);
basis = {'Taylor', 'Taylor'};
energy = 3;


%% set up the initial data

% source parameter in radians
localStableCoordinates = C.LocalMaps{1}(C.GlobalIntersection(1));
[s1_theta, s1_rad] = cart2pol(real(localStableCoordinates(1)), imag(localStableCoordinates(1)));  % retrieve arg(sigma) in the range [0, 2pi)
s1Init = s1_theta - pi; % translate arc(sigma) into [-pi, pi)
s1Init = s1Init - 0.3;
% s1_rad = s1_rad - 0.1;

% data for interior strands
cData = C.sample(1000);
u_data = cData.Orbit{1}(:, 1);  % u1 in X0
u_data(:, 2) = cData.Orbit{2}(:, 1);  % u2 in X1
u_data(:, 3) = cData.Orbit{3}(:, 1); % u3 in X0
u_data(:, 4) = cData.Orbit{4}(:, 1);  % u4 in X2
u_data(:, 5) = cData.Orbit{5}(:, 1);  % u5 in X0

% target parameter in radians
localUnstableCoordinates = C.LocalMaps{2}(C.GlobalIntersection(3));
[s2Theta, s2Rad] = cart2pol(real(localUnstableCoordinates(1)), imag(localUnstableCoordinates(1)));  % retrieve arg(sigma) in the range [0, 2pi)
s2Init = s2Theta - pi; % translate arg(sigma) into [-pi, pi)
s2Init = s2Init + 0.35  % correction for the bug in the local coordinate code

% fix all but 1 integration time and estimate the last one
fixedTau = arrayfun(@(j)sum(cData.Tau{j}), 1:4);
tauInit = sum(cData.Tau{5});

% initialize unfolding parameter at zero
alphaInit = 0;
uInit = [s1Init; reshape(u_data, [], 1); s2Init; tauInit; alphaInit];


%% set up the BVP map

% evaluate initial segment
Gu = [Scalar(Au, basis); Scalar(Bu, basis); Scalar(Cu, basis); Scalar(Du, basis)];
% Pu = cat(3, Au, Bu, Cu, Du);
Gamma_u = @(s)Gamma(Gu, s1_rad, s);
R1 = @(s, u1)(Gamma_u(s) - u1);


% vector fields and flows
f0 = @(t, u)rk45regvectorfield(t, u, mu, 0);  % physical vector field
f0_dissipative = @(t, u, alpha)(f0(t, u) + alpha*[0; u(2); 0; u(4)]);
f1 = @(t, u)rk45regvectorfield(t, u, [mu, energy], 1);  % regularized vector fields
f2 = @(t, u)rk45regvectorfield(t, u, [mu, energy], 2);
F0 = @(u, tau)tau_map(f0, [0, tau], u);
F1 = @(u, tau)tau_map(f1, [0, tau], u);
F2 = @(u, tau)tau_map(f2, [0, tau], u);
F0_alpha = @(u, tau, alpha)tau_map(@(t, u)f0_dissipative(t, u, alpha), [0, tau], u);

% conjugacy maps
h1Inverse = @(u)conj_inverse(u, mu, 1);
h2Inverse = @(u)conj_inverse(u, mu, 2);


% evaluate the interior shooting segments
R2 = @(u1, u2)h1Inverse(u2) - F0(u1, fixedTau(1));
R3 = @(u2, u3)h1Inverse(F1(u2, fixedTau(2))) - u3;
R4 = @(u3, u4)h2Inverse(u4) - F0(u3, fixedTau(3));
R5 = @(u4, u5)h2Inverse(F2(u4, fixedTau(4))) - u5;

% evaluate terminal segment
Gs = [Scalar(As, basis); Scalar(Bs, basis); Scalar(Cs, basis); Scalar(Ds, basis)];
% Ps = cat(3, As, Bs, Cs, Ds);
Gamma_s = @(s)Gamma(Gs, s2Rad, s);
R6 = @(u5, s, t, alpha)(F0_alpha(u5, t, alpha) - Gamma_s(s));
% u_data(:, 5) = Gamma_s(s2_init);

% vectorize the entire thing
F = @(u)[R1(u(1), u(2:5)); R2(u(2:5), u(6:9)); R3(u(6:9), u(10:13)); R4(u(10:13), u(14:17)); R5(u(14:17), u(18:21)); R6(u(18:21), u(22), u(23), u(24))];
DF = @(u)numeric_jacobian2(F, u, 1e-2);
% r = findroot(F, DF, u_init)

%% initialize Newton iteration
clc
x0 = uInit;
f = F;
df = DF;
% x = findroot(f, df, x0)

u1 = uInit(2:5);
u2 = uInit(6:9);
u3 = uInit(10:13);
u4 = uInit(14:17);
u5 = uInit(18:21);


y0 = f(x0); % initialize f(x)
Df = df(x0); % initialize Df(x)
x1 = x0 - Df\y0; % initialize Newton iteration

% return
%% plot orbit vs boundary of unstable vs point where they should meet
close all
figure
hold on
plot(cData.Orbit{1}(1,:), cData.Orbit{1}(3,:), 'b')
sFull = linspace(0, 2*pi);
unstable_bd = Gamma_u(sFull);
P1 = unstable_bd(1,:);
P3 = unstable_bd(3,:);
plot(P1, P3, 'b')
chk_local = Gamma_u(s1Init);
scatter(chk_local(1), chk_local(3))
plot_orbit(f0, [0, fixedTau(1)], cData.Orbit{1}(:, 1), [1,3])

% plot the same for the stable manifold
plot(cData.Orbit{end}(1,:), cData.Orbit{end}(3,:), 'r')
stable_bd = Gamma_s(sFull);
Q1 = stable_bd(1,:);
Q3 = stable_bd(3,:);
plot(Q1, Q3, 'r')
chk_local2 = Gamma_s(s2Init);
scatter(chk_local2(1), chk_local2(3))
plot_orbit(f0, [0, -tauInit], cData.Orbit{5}(:, end), [1,3])


% close all
figure
hold on
% plot(cData.Orbit{2}(1,:), cData.Orbit{2}(3,:), 'r')
% plot_orbit(f1, [0, fixedTau(2)], cData.Orbit{2}(:,1), [1,3], 'PlotOptions', {'g'});
% plot_orbit(f2, [0, fixedTau(2)], cData.Orbit{2}(:,1), [1,3], 'PlotOptions', {'b'});

atlas_ij = @(i,j)plot(cData.Orbit{2}(i,:), cData.Orbit{2}(j,:), 'r');
rk45_ij = @(i,j)plot_orbit(f1, [0, fixedTau(2)], cData.Orbit{2}(:,1), [i,j], 'PlotOptions', {'g'});
atlas_ij(1, 3)
rk45_ij(1, 3)


% for i = 1:3
%     for j = i+1:4
%         figure
%         hold on
%         atlas_ij(i, j)
%         rk45_ij(i, j)
%     end
% end

% compute trajectory using 0-dim taylor integrator
initialData = cData.Orbit{2}(:, 1);
% bd = RegCRTBPChart(initialData, 'Taylor', 0, 30, energy, [mu, energy], 1, 'InitialScaling', 0.1, 'boundary', true);
cht = RegCRTBPChart(initialData, 'Taylor', 0, 30, energy, [mu, energy], 1, 'InitialScaling', 0.1);
        
%         % Integrate multiple timesteps
%         bd = RegCRTBPChart(initialData, basis, initialTime, truncation, energy, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
%         maxTau = .3;
%         A1 = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
%         while ~isempty(A1.LeafStack)
%             A1.growboundary()
%         end
        
% plot orbits
% s = 0;

ob = cht.eval(linspace(0, 1, 100));
ob = cell2mat(ob')';
plot(ob(:,1), ob(:,3), 'k', 'LineWidth',2)

% set up rk45 evaluation
% VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
% initCondition = ob(1,:)';
% [~,sol] = ode45(VF,tf,initCondition,odeOptions);
% plot3(sol(:,1),sol(:,3),sol(:,2),'go')


dealfig()


function fout = Gamma(scalar, r, theta)
% write an evaluation function to wrap one of the local parameterizations

P = reshape(scalar, [], 1);  % coerce into a row of Scalars
z1 = r.*exp(1j*theta);
z2 = r.*exp(-1j*theta);
data = [z1, z2];
fout = cell2mat(P.eval(data));
end

function fout = conj_inverse(u, mu, conj_idx)
u = reshape(u, 1, []);  % CRTBP2reg requires a row vector
u0 = CRTBP2reg(u, mu, -conj_idx);
fout = reshape(u0(:, 1:4), [], 1); % chop off the automatic differentiation coordinates and turn into a column vector
end



