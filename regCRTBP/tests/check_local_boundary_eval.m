%CHECK_LOCAL_BOUNDARY_EVAL - Implement and test the local boundary for L4 using interpolation instead of polygonal decomposition

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 08-Mar-2023;

%% Get some data to test with
clearvars
close all
load newl4_equalmass_local_manifolds % get local manifold data

%%  first piece of the atlas run code

% ==================================== Analytic continuation parameters  ====================================
mu = 0.5; % mass ratio of small primary
C = 3; % Energy level set in which to search for connections
minTau = -2; % stable manifold integration time
maxTau = 2; % unstable manfiold integration time
N = 21; % spatial truncation
M = 40; % temporal truncation
truncation = [M, N]; % set truncation vector

% ==================================== Parameters to leave alone for the most part ====================================
regTime = @taylorregtime;  % set method of integration for regularizing time
subDivideParameter = [4, .1, 2]; % the new default is [4, .1, 2]
tauGuess = 0.1; % initial timestep to scale each Taylor step
[unstableBd, unstableLocalMap] = local_unstable_boundary(mu, C, 'newL4_equalmass_local_manifolds', 'newL4_equalmass_local_manifolds', truncation, tauGuess);
uBd = unstableBd.controltailratio(5, .1);
newBd = L4_boundary_via_chebyshev('newl4_equalmass_local_manifolds', truncation, 'unstable', C, tauGuess);

figure
hold on
sFull = linspace(-1, 1, 100).';


for j = 1:length(newBd)
    sc = newBd(j).Coordinate;
    Z = real(cell2mat(sc.eval(sFull)));
    scatter(Z(1,:), Z(3,:), 'filled')
end

for j = 1:length(uBd)
    sc = uBd(j).Coordinate;
    Z = real(cell2mat(sc.eval(sFull)));
    plot(Z(1,:), Z(3,:), 'k', 'LineWidth', 2)
end