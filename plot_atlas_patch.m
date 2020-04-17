%PLOT_ATLAS - script for testing atlas plotting and temporarily for
%plotting atlases as patches.
%
%   Description:
%       PLOT_ATLAS description
%
%   Output:
%       PLOT_ATLAS output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 20-Sep-2019; Last revision: 20-Sep-2019

%% ================================================== PLOT AN ATLAS ==================================================
close all
% this assumes there is an object called plotAtlas in the workspace
if ~exist('plotIdx')
    plotIdx = 1:plotAtlas.Size;
end

% set up plot
figure;
hold on

% set up an evaluation grid
tValue = linspace(0,tau,50)';
tGlobal = linspace(0.1, plotAtlas.MaxTau, 200)';
sGlobal = linspace(-1, 1, 25)';
[S,T] = meshgrid(sGlobal, tGlobal);
evalData = [reshape(S, [],1), reshape(T, [], 1)];

% compute patch data
face = [];
vertex = [];
colorData = [];
for jChart = plotAtlas.Chart(plotIdx)
    bdVertices = [jChart.SpatialSpan(1), jChart.TimeSpan(1); jChart.SpatialSpan(1), jChart.TimeSpan(2);... % add 4 corners of domain
        jChart.SpatialSpan(2), jChart.TimeSpan(1); jChart.SpatialSpan(2), jChart.TimeSpan(2)];
    data = [jChart.intersectdomain(evalData); bdVertices];
    evalChart = cell2mat(jChart.eval(data, 'GlobalTime', true, 'GlobalSpace', true));
    
    % map evaluations to F0 coordinates
    if ~isequal(jChart.RegType, 0)
        mu = jChart.Parameter(1);
        F0Data = CRTBP2reg(evalChart(:,1:4), mu, -jChart.RegType);
    else
        F0Data = evalChart;
    end
    
    % choose coordinates to plot
    x = F0Data(:,1);
    y = F0Data(:,3);
    z = F0Data(:,2);
        
    jColorData = data(:,2);
    jVertex = [mid(x), mid(y) ,mid(z)];
    jFace = delaunay(data(:,1), data(:,2)); % index triples for delaunay triangulation in space-time.
    vertexCount = size(vertex,1);
    face = [face; jFace + vertexCount]; % start counting face labels from the previous end label
    vertex = [vertex; jVertex]; % append new vertices
    if ~isequal(length(jColorData), size(jVertex,1))
        disp('something went wrong here')
    end
    colorData = [colorData;jColorData];
end

% plot the patch
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,  'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha',1);
colorbar()
dealfig() 

% plot bounding circle
% bdtheta = linspace(0,2*pi,500);
% xb = 0.5*cos(bdtheta);
% yb = 0.5*sin(bdtheta);
% pb = zeros(size(bdtheta));
% plot3(xb,yb,pb,'y','LineWidth',3)


% view(10,9)
% set(gca,'ZLim',[-15,15])



return
%% ======================= CODE SNIPPET FOR PLOTTING ATLASES WITH SCATTER =====================

% plot interior and initial boundary
for jChart = plotAtlas.Chart(plotIdx)
    evalInterior = cell2mat(jChart.eval(evalData));
    
    % choose coordinates to plot
    x = evalInterior(:,1);
    y = evalInterior(:,3);
    z = evalInterior(:,2);
    
    if jChart.Generation == 0
        evalBoundary = cell2mat(jChart.eval([sLocal, zeros(size(sLocal))]));
        x0 = evalBoundary(:,1);
        y0 = evalBoundary(:,3);
        z0 = evalBoundary(:,2);
        plot3(mid(x0), mid(y0), mid(z0), 'k', 'LineWidth',3)
    end
    scatter3(mid(x), mid(y), mid(z), 5,'filled', 'b')
end


% plot terminal boundaries
% termChart = RegCRTBPChart;
% for jChart = 1:length(plotAtlas.CrashStack)
%     parentChart = plotAtlas.CrashStack(jChart).ParentHandle;
%     termChart(jChart) = parentChart;
%     [x1,p1,y1,q1,z1] = parentChart.eval([sRelative, ones(size(sRelative))]);
%     plot3(mid(x1),mid(y1),mid(p1), 'r', 'LineWidth',3)
% end

% plot orbits
% tf = linspace(0, plotAtlas.MaxTau, 100);
% for s = [-1,1]
%     ob = mid(plotAtlas.orbit(s, tf));
%     plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',3)
% end

% tb = linspace(0,backwardTau,100);
% ob = mid(backwardAtlas.orbit(0, tb));
% plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',3)




