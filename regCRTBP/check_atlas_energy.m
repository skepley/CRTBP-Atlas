%CHECK_ATLAS_ENERGY - Compute energy statistics for all charts in an atlas.
%
%   Description:
%       CHECK_ATLAS_ENERGY description
%
%   Output:
%       CHECK_ATLAS_ENERGY output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 08-Jul-2019; Last revision: 08-Jul-2019

% load([savePath, 'P1_equalmass_L4energy_t5backward'])
% chkAtlas = atlasP1B;
% load([savePath, 'L4_equalmass_t5backward'])
chkAtlas = atlasL4B;

% get statistical data
nPts = 10;
[chartEnergy, chartSTD] = arrayfun(@(chart)chart.energystats(nPts), chkAtlas.Chart);
regData = [chkAtlas.Chart.RegType];
timeData = [chkAtlas.Chart.TimeSpan];
midTimeData = (timeData(1:2:end) + timeData(2:2:end))/2;

meanEnergy = mean(chartEnergy); % mean energy over all charts of this atlas
stdEnergy = std(chartEnergy); % standard deviation of energy over all charts of this atlas
maxEnergy = max(chartEnergy);
minEnergy = min(chartEnergy);

% % identify bad charts and call plot_atlas_patch
badChartIdx = abs(chartEnergy-meanEnergy) > 3*stdEnergy;
plotAtlas = chkAtlas;
% plotIdx = badChartIdx;
plotIdx = 1:chkAtlas.Size;

%% ================================================== PLOT THE ATLAS ==================================================
% plot the entire manifold with color showing deviation from mean energy

% set up plot
% close all
figure;
hold on

% plot primaries
primarySize = 500;
scatter3(mu, 0, 0, (1-mu)*primarySize,'filled','b')
scatter3(mu-1, 0 , 0, mu*primarySize,'filled','k')
scatter3(L4(1), L4(3),0, primarySize, 'rx')

% plot bounding circles
bdtheta = linspace(0,2*pi,500);
xb = 0.25*cos(bdtheta);
yb = 0.25*sin(bdtheta);
pb = zeros(size(bdtheta));
plot3(mu + xb,yb,pb,'r','LineWidth',2)
plot3(mu-1 + xb,yb,pb,'r','LineWidth',2)

% set up an evaluation grid
tGlobal = linspace(0.1, chkAtlas.MaxTau, 500)';
sGlobal = linspace(-1, 1, 250)';
[S,T] = meshgrid(sGlobal, tGlobal);
evalData = [reshape(S, [],1), reshape(T, [], 1)];

% compute patch data
face = [];
vertex = [];
colorData = [];
for jChart = chkAtlas.Chart(plotIdx)
    bdVertices = [jChart.SpatialSpan(1), jChart.TimeSpan(1); jChart.SpatialSpan(1), jChart.TimeSpan(2);... % add 4 corners of domain
        jChart.SpatialSpan(2), jChart.TimeSpan(1); jChart.SpatialSpan(2), jChart.TimeSpan(2)];
    data = [jChart.intersectdomain(evalData); bdVertices];
    %     data = jChart.intersectdomain(evalData);
    evalChart = cell2mat(jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true));
    
    % map evaluations to F0 coordinates
    if ~isequal(jChart.RegType, 0)
        mu = jChart.Parameter(1);
        F0Data = CRTBP2reg(evalChart, mu, -jChart.RegType);
    else
        F0Data = evalChart;
    end
    
    if ~isempty(F0Data)
        % choose coordinates to plot
        x = F0Data(:,1);
        y = F0Data(:,3);
        z = F0Data(:,2);
        
        % patch color by energy
        jEnergy = jChart.energystats(nPts);
        jColorData = min(jEnergy, 3.5)*ones(size(x));
        
        % patch color by ell^1 norm
        %     jColorData = min(max(jChart.ellonebox.rad), 1)*ones(size(x));
        try
            jVertex = [mid(x), mid(y) ,mid(z)];
            jFace = delaunay(data(:,1), data(:,2)); % index triples for delaunay triangulation in space-time.
            vertexCount = size(vertex,1);
            face = [face; jFace + vertexCount]; % start counting face labels from the previous end label
            vertex = [vertex; jVertex]; % append new vertices
            if ~isequal(length(jColorData), size(jVertex,1))
                disp('something went wrong here')
            end
            colorData = [colorData;jColorData];
        catch
            
        end
    end
end

% plot the patch
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
title(sprintf('sd: %0.3f', subDivideParameter(2)));
colorbar()
dealfig()



%% ================================================== CHECK WITH RK45 ==================================================

regType = 0;
F = @(t,x)rk45regvectorfield(t,x,mu,regType);

% set up an evaluation grid
% tLocal = linspace(0, 1, 10)';
sLocal = linspace(-1, 1, 3)';
% [S,T] = meshgrid(sLocal, tLocal);
% evalData = [reshape(S, [],1), reshape(T, [], 1)];
evalData = [sLocal, zeros(size(sLocal))];
genZeroCount = 0;
odeOptions = odeset('RelTol',1e-14,'AbsTol',1e-15);
for jChart = atlasL4B.Chart
    if isequal(jChart.Generation, 0)
        genZeroCount = genZeroCount + 1;
        evalChart = cell2mat(jChart.eval(evalData));
        for j = 1:size(evalChart,1)
            ic = evalChart(j,:)';
            [~,jOrbit] = ode45(F, [0,maxTau], ic, odeOptions);
            plot3(jOrbit(:,1), jOrbit(:,3), jOrbit(:,2), 'r')
            jEnergy = CRTBPenergy(jOrbit', mu, 0);
            disp(max(jEnergy))
%             plot_orbit(F, [0,maxTau], ic, [1,3,2], 'PlotOptions', {'r'})
        end
    end
end

%% ===========================PLOT ENERGY STATS FOR ATLAS ===============================
return

figure
scatter(midTimeData, chartSTD, 5, 'filled')
title('std')

figure
scatter(midTimeData, chartEnergy, 5, 'filled')
title('mean')

maxNorm = zeros(size(midTimeData));
allNorm = zeros(6, length(midTimeData));
chartCenter = zeros(length(midTimeData), 2);
for j = 1:chkAtlas.Size
    allNorm(:,j) = chkAtlas.Chart(j).norm;
    maxNorm(j) = max(allNorm(:,j));
    chartCenter(j, :) = [chkAtlas.Chart(j).Coordinate(1).Coefficient(1), chkAtlas.Chart(j).Coordinate(3).Coefficient(1)];
end

figure
scatter(midTimeData, maxNorm, 5, 'filled')
title('max norm')

figure
hold on
for kk = 1:4
    scatter(midTimeData, allNorm(kk,:), 5, 'filled')
end
legend()

figure
hold on
scatter(chartCenter(:,1), chartCenter(:,2), 5, 'filled')

% plot primaries
primarySize = 500;
scatter(mu,0,(1-mu)*primarySize,'filled','b')
scatter(mu-1,0,mu*primarySize,'filled','k')
scatter(L4(1), L4(3), primarySize, 'rx')

% figure
% scatter(1:length(regData), regData, 5, 'filled')
% title('reg')


dealfig()

