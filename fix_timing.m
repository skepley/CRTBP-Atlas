%FIX_TIMING - Figure out the problem with the timing. This seems to occur when a chart is in regtype 1 or 2 coordinates
%and reaches the terminal time for the computation. At this point the time rescaling is doing something wonky and the
%final chart has a norm ~1e11
%
%   Description:
%       FIX_TIMING description
%
%   Output:
%       FIX_TIMING output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 19-May-2021; Last revision: 19-May-2021

%% ================================================== SECTION 1 ==================================================
% Start by loading the 7 time unit L4 unstable manifold as atlasL4F
clear C
clc
atl = atlasL4F;
B = atl.Chart(14142); % this is one of the bad charts
A = B.ParentHandle; % its parent chart is fine.

% Lets walk through the construction of B from A defined as C
C = A.fixtime(1);
tau = atl.TimeStepper(C); % get tau from timestepper
C.advect(tau); % convert boundary chart to interior chart by advection
% disp(C);
C.rescaletime(eps(1)); % rescale timestep to force last coefficient below a given norm
% disp(C);
taylorregtime(C); % normalize time to the global F0 time scale
% disp(C);


% next round
D = C.fixtime(1);  % new boundary
tau = atl.TimeStepper(D); % get tau from timestepper
D.advect(tau); % convert boundary chart to interior chart by advection
D.rescaletime(eps(1)); % rescale timestep to force last coefficient below a given norm
disp(D);
taylorregtime(D); % normalize time to the global F0 time scale

% now the block from boundarycheck
maxTau = 7;
newTau = maxTau - D.TimeSpan(1); % get timestep in f0 time
newRegTime = inverseregtime(D, newTau); % compute corresponding timestep in F1/F2 time
disp(D)
D.Tau = D.LocalTime;
D.scaletime(newRegTime);  % rescale time so that max time is exactly maxTau
a = rk45regtime(D)
b = taylorregtime(D, D.Tau)
taylorregtime(D);
% disp(D);


%% Use taylor to get F0 time which should be exactly equal to newTau
F1Tau = D.Tau;
x_i = D.Coordinate(1).Coefficient(:,1); % Taylor expansion in time for the orbit through x_i(0)
y_i = D.Coordinate(3).Coefficient(:,1); % Taylor expansion in time for the orbit through y_i(0)
M = D.Truncation(1); % Truncation size in time
% t = 4*tau*(sum(conv(x_i, x_i)./(1:2*M-1).') + sum(conv(y_i, y_i)./(1:2*M-1).')); % integrate dt/dTau_i

p = fliplr(conv(x_i, x_i).');
q = fliplr(conv(y_i, y_i).');
t = 4*F1Tau*(polyval(polyint(p + q), D.Tau/F1Tau));
disp(t - newTau);

% use rk45 to get F0 time which is exactly .0026 as it should be
initialCondition = cell2mat(D.eval(zeros(1, D.Dimension(1)))); % get an initial point on the chart
F = @(t,x)rk45regvectorfield(t, x, D.Parameter, D.RegType);
[tau, Fsol] = ode45(F, D.TimeSpan, initialCondition.', odeset('AbsTol', 1e-13));
w = Fsol(:,1).^2 + Fsol(:,3).^2;
F0Tau = 4*trapz(tau,w);
disp(F0Tau - newTau);

F0Tau2 = integral(@(tau_vector)arrayfun(@(tau)g(F, initialCondition, [D.TimeSpan(1), tau]), tau_vector), D.TimeSpan(1), D.TimeSpan(2));
disp(F0Tau2 - newTau);



%% So the Taylor orbits are somehow wrong? 
nNode = 100;
Xt = polyval(fliplr(x_i.'), linspace(0, 1, nNode));
Yt = polyval(fliplr(y_i.'), linspace(0, 1, nNode));
rkTimeNodes = linspace(0, D.LocalTime, nNode);
[tau, Fsol] = ode45(F, rkTimeNodes, initialCondition.', odeset('AbsTol', 1e-13));
Xr = Fsol(:, 1);
Yr = Fsol(:, 3);

close all
figure
hold on
plot(Xt, Yt, 'r', 'LineWidth', 3);
plot(Xr, Yr, 'g', 'LineWidth', 1);





function fout = g(F, IC, f1_time)
% a callable function to use in the MATLAB integral routine
    xf = tau_map(F, f1_time, IC.');
    fout = 4*(xf(1).^2 + xf(3).^2);
end
% return
% % now it crashed
% disp(D);
% E = D.fixtime(1);
% boundarycheck(atl, E, 7)
% disp(D)
% disp(E)


