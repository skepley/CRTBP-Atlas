%CHECK_CONNECTION

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 10-Aug-2020; Last revision: 19-Aug-2020


% first run final2.m
close all
% clearvars
computerPath = pcpath('mac');
addpath('/users/shane/dropbox/matlab/tools/exportfig')
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
load newl4_equalmass_local_manifolds % get local manifold data
load chk_connection_save.mat
% clear connections
% for j=length(sols):-1:1
%     connections(j) = RegCRTBPConnection(sols{j}, stableLocalMap, unstableLocalMap);
% end
% ctimes = [connections.ConnectionTime];
% min(ctimes)
% save('chk_connection_save.mat', 'connections')
mu = connections(1).Parameter;

return
%% test functionality of RegCRTBPConnection class


% j = 29;
% C = RegCRTBPConnection(sols{j}, stableLocalMap, unstableLocalMap);
% figure('OuterPosition', [10, 10, 1600, 900])
% hold on
% C.plot(1000)
% C.plot_connection_tangents()
% rk45_from_intersection(C)



%% Try to verify connections using RK45
close all
trueSolIdx = find([connections.TrueOrbit] == true);
% plotIdx = 1:10;  % look at 20 orbits at a time.

figure 
hold on
for k = trueSolIdx
    disp(k)
%     figure('OuterPosition', [10, 10, 1600, 900])
%     hold on
    connections(k).plot(1000, 'LineWidth', 1.5)
    %     connections(k).plot_connection_tangents()
    %     try
    %         rk45_from_intersection(connections(k))
    %         plot_primaries(mu, L4)
    %         title(sprintf('%d Score: %0.8f', k, rk45_score(connections(k))))
    %         axis equal
    %         export_fig('-png', '-transparent', sprintf('/users/shane/dropbox/Regularisation3bp/CRTBP Atlas/production runs/connections/%d.png', k))
    %         close all
    %     catch
    %         warning(sprintf('Caught one at %d',k))
    %     end
    
end
return
% %% Try to verify connections using tangency check at the intersection point
% % clear connections
% % for j=length(sols):-1:1
% %     connections(j) = RegCRTBPConnection(sols{j}, stableLocalMap, unstableLocalMap);
% % end
%
% close all
% trueSolIdx = find([connections.TrueOrbit] == true);
% plotIdx = 1:10;  % look at 20 orbits at a time.
%
% tangentSolCheck = @(C)(tangent_score(C) > 0.98);
%
% for k = 1:length(connections)
%     if isequal(connections(k).StableChart.RegType, connections(k).UnstableChart.RegType) && tangentSolCheck(connections(k))
%         disp(tangentSolCheck(connections(k)))
%
%         figure('OuterPosition', [10 + 2000, 10, 1600, 900])
%         hold on
%         connections(k).plot(10000, 'LineWidth', 1.5)
%         connections(k).plot_connection_tangents()
%         plot_primaries(mu, L4)
%         title(sprintf('%d Score: %0.8f', k, tangent_score(connections(k))))
%         axis equal
%         export_fig('-png', sprintf('/users/shane/dropbox/Regularisation3bp/CRTBP Atlas/production runs/tangent_connections/%d.png', k))
%         close all
%     end
%
% end

%% These are some pretty ones
prettyOnes = [45, 193, 461,464, 469];

for k = prettyOnes
    figure('OuterPosition', [10, 10, 1600, 900])
    hold on
    connections(k).plot(50000, 'LineWidth', 1.5)
    plot_primaries(mu, L4)
    axis equal
    export_fig('-png', sprintf('/users/shane/dropbox/Regularisation3bp/CRTBP Atlas/production runs/connections/pretty_%d.png', k))
    dealfig()
end


%% CODE HUNK FOR CHECKING ORBITS AGAINST RK45
%     F0 = @(t,x)rk45regvectorfield(t, x, mu, 0);  % rk45 field to test against

%     F0InitialData = C.Orbit{1}(:,1);
%     F0FinalData = C.Orbit{end}(:,end);

%     r0 = 1./sqrt((F0InitialData(1) - mu).^2 + F0InitialData(3).^2);
%     s0 = 1./sqrt((F0InitialData(1) + mu).^2 + F0InitialData(3).^2);
%     F0InitialData = [F0InitialData; r0; s0];

%     plot_orbit(F0, [0, C.UnstableTime], addAutoDiff(F0InitialData, mu), [1,3], 'PlotOptions', {'g-.', 'LineWidth', 2})
%     plot_orbit(F0, [0, C.UnstableTime], addAutoDiff(F0InitialData, mu), [1,3], 'PlotOptions', {'k-.', 'LineWidth', 3})
%     plot_orbit(F0, [0, -(C.ConnectionTime - C.UnstableTime)], addAutoDiff(F0FinalData, mu), [1,3],'PlotOptions', {'k-.', 'LineWidth', 2})



function plot_primaries(mu, L4)
% plot primaries and L4
primarySize = 100;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','k')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')
plot3(L4(1), L4(3), L4(2), 'r*')
end

% function score = tangent_score(C)
% % score connection based on how close to tangent the vector field is at the two intersection images
% mu = C.Parameter;
% if isequal(C.StableChart.RegType, C.UnstableChart.RegType)
%     regType = C.StableChart.RegType;
% else
%     error('This is not implemented yet')
% end
%
% F = @(t,x)rk45regvectorfield(t, x, mu, regType);  % rk45 field to test against
% stableInitialData = cell2mat(C.StableChart.eval(C.LocalIntersection(1:2).'));
% unstableInitialData = cell2mat(C.UnstableChart.eval(C.LocalIntersection(3:4).'));
% Fs = F(0, stableInitialData);
% Fu = F(0, unstableInitialData);
% score = dot(Fs, Fu)/(norm(Fs)*norm(Fu));  % return cos(theta) where theta is the angle between Fu and Fs
% end

function rk45_from_intersection(C)
mu = C.Parameter;
if isequal(C.StableChart.RegType, C.UnstableChart.RegType)
    regType = C.StableChart.RegType;
else
    error('This is not implemented yet')
end

F = @(t,x)rk45regvectorfield(t, x, mu, regType);  % rk45 field to test against
stableInitialData = cell2mat(C.StableChart.eval(C.LocalIntersection(1:2).'));
unstableInitialData = cell2mat(C.UnstableChart.eval(C.LocalIntersection(3:4).'));

plot_orbit(F, [0, 0.5], stableInitialData, [1,3], 'PlotOptions', {'b.', 'LineWidth', 2})
plot_orbit(F, [0, 0.5], unstableInitialData, [1,3], 'PlotOptions', {'k.', 'LineWidth', 2})
plot_orbit(F, [0, -0.5], stableInitialData, [1,3], 'PlotOptions', {'b.', 'LineWidth', 2})
plot_orbit(F, [0, -0.5], unstableInitialData, [1,3], 'PlotOptions', {'k.', 'LineWidth', 2})
end

function score = rk45_score(C)
% score connection based on far the orbits are after +- 0.5 time units after integrating from the two intersection images

mu = C.Parameter;
if isequal(C.StableChart.RegType, C.UnstableChart.RegType)
    regType = C.StableChart.RegType;
else
    error('This is not implemented yet')
end
F = @(t,x)rk45regvectorfield(t, x, mu, regType);  % rk45 field to test against
stableInitialData = cell2mat(C.StableChart.eval(C.LocalIntersection(1:2).'));
unstableInitialData = cell2mat(C.UnstableChart.eval(C.LocalIntersection(3:4).'));
forward_score = norm(tau_map(F, [0, 0.5], stableInitialData) - tau_map(F, [0, 0.5], unstableInitialData));
backward_score = norm(tau_map(F, [0, -0.5], stableInitialData) - tau_map(F, [0, -0.5], unstableInitialData));
score = forward_score + backward_score;
end

% function uFull = addAutoDiff(u, mu)
% r0 = 1./sqrt((u(1) - mu).^2 + u(3).^2);
% s0 = 1./sqrt((u(1) + mu).^2 + u(3).^2);
% uFull = [u; r0; s0];
% end

