%CHECK_CONNECTION

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 10-Aug-2020; Last revision: 19-Aug-2020
load newl4_equalmass_local_manifolds % get local manifold data



clear connections
nNode = 1000;
trueSols = [];
for j=1:length(sols)
    if exist('connections')
        try
            connections(end+1) = connection2BVP(sols{j}, nNode, stableLocalMap, unstableLocalMap);
            trueSols(end+1) = sols{j}{4};
        catch
            sprintf('Bad sol = %d', j)
        end
        
    else
        try
            connections = connection2BVP(sols{j}, nNode, stableLocalMap, unstableLocalMap);
            trueSols(end+1) = sols{j}{4};
            
        catch
            sprintf('Bad sol = %d', j)
        end
        
    end
end


%%

close all
F0 = @(t,x)rk45regvectorfield(t, x, mu, 0);  % rk45 field to test against
trueSolIdx = find(trueSols == true);

figure
hold on
% all orbits in 1 plot
for k = trueSolIdx
    disp(k)
    C = connections(k);
    disp(C)
%     figure
    hold on
    for j = 1:length(C.Orbit)
        switch C.RegVector(j)
            case 0
                jSeg = C.Orbit{j};
                plot(jSeg(1,:), jSeg(3,:), 'g', 'LineWidth', 1)

            case 1
                jSeg = CRTBP2reg(C.Orbit{j}.', mu, -1).';
                plot(jSeg(1,:), jSeg(3,:), 'r', 'LineWidth', 1)

            case 2
                jSeg = CRTBP2reg(C.Orbit{j}.', mu, -2).';
                plot(jSeg(1,:), jSeg(3,:), 'b', 'LineWidth', 1)
        end
    end
    F0InitialData = C.Orbit{1}(:,1);
    F0FinalData = C.Orbit{end}(:,end);

    %     r0 = 1./sqrt((F0InitialData(1) - mu).^2 + F0InitialData(3).^2);
    %     s0 = 1./sqrt((F0InitialData(1) + mu).^2 + F0InitialData(3).^2);
    %     F0InitialData = [F0InitialData; r0; s0];

    %     plot_orbit(F0, [0, C.UnstableTime], addAutoDiff(F0InitialData, mu), [1,3], 'PlotOptions', {'g-.', 'LineWidth', 2})
    %     plot_orbit(F0, [0, C.UnstableTime], addAutoDiff(F0InitialData, mu), [1,3], 'PlotOptions', {'k-.', 'LineWidth', 3})
    %     plot_orbit(F0, [0, -(C.ConnectionTime - C.UnstableTime)], addAutoDiff(F0FinalData, mu), [1,3],'PlotOptions', {'k-.', 'LineWidth', 2})
end

% return
% % separate plots
% for k = trueSolIdx
%     C = connections(k);
%     pl_k = figure;
%     hold on
%     for j = 1:length(C.Orbit)
%         switch C.RegVector(j)
%             case 0
%                 jSeg = C.Orbit{j};
%                 plot(jSeg(1,:), jSeg(3,:), 'g', 'LineWidth', 1)
%                 
%             case 1
%                 jSeg = CRTBP2reg(C.Orbit{j}.', mu, -1).';
%                 plot(jSeg(1,:), jSeg(3,:), 'r', 'LineWidth', 1)
%                 
%             case 2
%                 jSeg = CRTBP2reg(C.Orbit{j}.', mu, -2).';
%                 plot(jSeg(1,:), jSeg(3,:), 'b', 'LineWidth', 1)
%         end
%     end
%     plot_primaries(mu, L4)
%     title(sprintf('%d', k))
%     saveas(pl_k, sprintf('/users/shane/desktop/homoclinics/orbit_%d.png', k))
%     close 1
% end


%
trueConnections = 5;
plot_primaries(mu, L4)




function vecRep = connection_vector(con)
vecRep = [con.ConnectionTime, localUnstable(1), localStable(1)];  % representation to compare connection properties as vectors in C^4.
end

function plot_primaries(mu, L4)
% plot primaries and L4
primarySize = 100;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','k')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')
plot3(L4(1), L4(3), L4(2), 'r*')
end

function uFull = addAutoDiff(u, mu)
r0 = 1./sqrt((u(1) - mu).^2 + u(3).^2);
s0 = 1./sqrt((u(1) + mu).^2 + u(3).^2);
uFull = [u; r0; s0];
end


