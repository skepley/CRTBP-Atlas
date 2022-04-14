%CHECK_CONNECTION_CLASS - Import a file of connection data and extract and plot the oribts using rk45 and the Taylor charts
%
%   Description:
%       CHECK_CONNECTION_CLASS description
%
%   Output:
%       CHECK_CONNECTION_CLASS output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jun-2019; Last revision: 18-Jun-2019

%% process solutions

% clear all
% close all
% clc
% computerPath = pcpath('mac');
% addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
% addpath(genpath('/Users/sk2011/Dropbox/Regularisation3bp/integrator'))
% savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];
% 


% %% LOAD P1 TO P2 HETEROCLINICS
% load([savePath, 'P1_equalmass_L4energy_t5backward'])
% atlasB = atlasP1B;
% load([savePath, 'P2_equalmass_L4energy_t5forward'])
% atlasF = atlasP2F;
% load([savePath, 'P1P2_collision_heteroclinics_equalmass_L4energy'])
% 
% %% LOAD P1 TO P1 HOMOCLINICS
% load([savePath, 'P1_equalmass_L4energy_t5backward'])
% atlasB = atlasP1B;
% load([savePath, 'P1_equalmass_L4energy_t5forward'])
% atlasF = atlasP1F;
% load([savePath, 'P1_collision_homoclinics_equalmass_L4energy'])
% 
% %% LOAD L4 HOMOCLINICS
% load([savePath, 'L4_equalmass_t5backward'])
% atlasB = atlasL4B;
% load([savePath, 'L4_equalmass_t5forward'])
% atlasF = atlasL4F;
% load([savePath, 'L4_homoclinics_equalmass'])

%%
% close all
% figure
% hold on

connections = [sols{:}];
for jConn = 1:1 % size(connections,2)
    disp(jConn)
    cData = connections(:,jConn);
    
    chartB = cData{1}; % backward time chart
    chartF = cData{2}; % forward time chart
    localIntersection = cat(1, cData{3}(1:2).',  [cData{3}(3), 0]);
    globalIntersection = cat(1, chartB.local2global(localIntersection(1,:)), chartF.local2global(localIntersection(2,:)));
    phaseIntersection = cell2mat(chartB.eval(localIntersection(1,:))); % A point in R^n which lies on the connection
    
    % swap to F0 field
    switch chartB.RegType
        case 0 % do nothing
        case 1
            phaseIntersection = CRTBP2reg(phaseIntersection(1:4), chartB.Parameter(1), -1);
        case 2
            phaseIntersection = CRTBP2reg(phaseIntersection(1:4), chartB.Parameter(1), -2);
    end % switch
    
    phasePoint = phaseIntersection(1:6);
    flightTime = abs(diff(globalIntersection(:,2)));
    nPoint = 5000; % number of sample points on each orbit
    
    flightTimeB = min(globalIntersection(:,2));
    flightTimeF = max(globalIntersection(:,2));
    
    minTime = 1e-3; % for t = 0 both orbits are undefined in F0 space so set this minimum > 0.
    nB = round(-flightTimeB*nPoint/flightTime); % number of samples in backward time
    tB = linspace(-minTime, flightTimeB, nB)'; % backward time samples
    nF = nPoint - nB; % number of samples in forward time
    tF = linspace(minTime, flightTimeF, nF)';
    ob1 = mid(CRTBPorbit(atlasB, globalIntersection(1,1), tB));
    ob2 = mid(CRTBPorbit(atlasF, globalIntersection(2,1), tF));
    
    % reorient into increasing order
    ob = cat(1, flip(ob1,1), ob2);
    t = cat(1, flip(tB,1), tF);
    
    %  plot connections
    lineWidth = 1;
    plot(ob1(:,1), ob1(:,3), 'Color', [1,0,0], 'LineWidth', lineWidth + 1)
%     scatter(ob1(:,1), ob1(:,3))
    plot(ob2(:,1), ob2(:,3), 'Color', [1,0,0], 'LineWidth', lineWidth+1)
    %     scatter(ob2(:,1), ob2(:,3))
    
%     scatter(phasePoint(1), phasePoint(3), 'b*')

    
    
    
    % PLOT USING ODE 45
    %     F0 = @(t,x)rk45regvectorfield(t,x,mu,0);
    %     plot_orbit(@(t,x)-F0(t,x), [0, -.9*globalIntersection(1,2)], phasePoint, [1,3], 'PlotOptions', {'Color', [1,0,0], 'LineWidth', lineWidth })
    %     plot_orbit(F0, [0, .9*globalIntersection(2,2)], phasePoint, [1,3], 'PlotOptions', {'Color', [0,1,0], 'LineWidth', lineWidth })
    
end

% plot primaries
primarySize = 500;
scatter(mu,0,(1-mu)*primarySize,'filled','b')
scatter(mu-1,0,mu*primarySize,'filled','k')
scatter(L4(1), L4(3), primarySize, 'rx')

dealfig()
return


%% PLOT L4 LOCAL MANIFOLDS

% plot the manifold in F0 coordinates
s = linspace(-1,1,10);
tFinal = atlasL4F.MaxTau;
t = linspace(0, 1, 10);
[S,T] = meshgrid(s,t);
evalData = [reshape(S,[],1), reshape(T,[],1)];
close all
figure
hold on

for jChart = atlasL4F.Chart
    
    % map evaluations to F0 coordinates
    if isequal(jChart.RegType, 1)
        [X1,P1,Y1,Q1,~] = jChart.eval(evalData);
        F0Data = CRTBP2reg([X1,P1,Y1,Q1],mu,-1);
    elseif isequal(jChart.RegType, 2)
        [X2,P2,Y2,Q2,~] = jChart.eval(evalData);
        F0Data = CRTBP2reg([X2,P2,Y2,Q2],mu,-2);
    else
        [X0, P0, Y0, Q0, ~, ~] = jChart.eval(evalData);
        F0Data = [X0, P0, Y0, Q0];
    end
    
    X = reshape(F0Data(:,1), length(s), []);
    Y = reshape(F0Data(:,3), length(s), []);
    Z = reshape(F0Data(:,2), length(s), []);
    surf(X,Y,Z)
    shading interp
end
colorbar






%% plot local manifold interior
interior = CMTimestep();
interior.Coord = BAscalar(squeeze(P(1,:,:)));
for j = 1:4
    interior.Coord(j) = BAscalar(squeeze(P(j,:,:)));
end
u = interior.eval(interiorNode); % cooresponding coordinates for each vertex

%% FIX EVALUATION AND PLOTTING BUG
cB = cData{1};
cF = cData{2};
sts = cData{3};
gI = globalIntersection;
lI = localIntersection;
sB = lI(1,1);
tB = lI(1,2);

% eval 1
s = linspace(-1, sB, 50);
t = linspace(0,tB,50);
[S,T] = meshgrid(s,t);
data = [reshape(S, [], 1), reshape(T, [], 1)];
evl1 = cB.eval(data);

% eval 2
s = linspace(sB, 1, 50);
t = linspace(0,tB,50);
[S,T] = meshgrid(s,t);
data = [reshape(S, [], 1), reshape(T, [], 1)];
evl2 = cB.eval(data);

% eval 3
s = linspace(-1, sB, 50);
t = linspace(tB,1, 50);
[S,T] = meshgrid(s,t);
data = [reshape(S, [], 1), reshape(T, [], 1)];
evl3 = cB.eval(data);

% eval 4
s = linspace(sB, 1, 50);
t = linspace(tB, 1,50);
[S,T] = meshgrid(s,t);
data = [reshape(S, [], 1), reshape(T, [], 1)];
evl4 = cB.eval(data);

% orbit segment
t = linspace(cB.TimeSpan(1), gI(1,2), 10);
obSeg = mid(CRTBPorbit(atlasB, gI(1,1), t));

figure
hold on
scatter(evl1{1}, evl1{3}, 5, 'filled')
scatter(evl2{1}, evl2{3}, 5, 'filled')
% scatter(evl3{1}, evl3{3}, 5, 'filled')
% scatter(evl4{1}, evl4{3}, 5, 'filled')
plot(obSeg(:,1), obSeg(:,3), 'Color', [0,1,0], 'LineWidth', lineWidth)
scatter(phasePoint(1), phasePoint(3), 'b*')
dealfig()




return


%% FIND DUPLICATES
% close all
clc
C1 = connections(1,:); C1 = [C1{:}];
C2 = connections(2,:); C2 = [C2{:}];
ID = @(chart,j)chart.InitialData(j).Coefficient;
IDn = @(chart1, chart2, j)norm(ID(chart1, j) - ID(chart2,j),1);
IDv = @(chart1,chart2)norm(arrayfun(@(j)IDn(chart1,chart2,j),1:chart1.Dimension(2)),1);

A1 = C1(2);
B1 = C1(3);
IDv(A1,B1)

A2 = A1.ParentHandle;
B2 = B1.ParentHandle;
IDv(A2,B2)

A3 = A2.ParentHandle;
B3 = B2.ParentHandle;
IDv(A3,B3)

A4 = A3.ParentHandle;
B4 = B3.ParentHandle;
IDv(A4,B4)

A5 = A4.ParentHandle;
B5 = B4.ParentHandle;
IDv(A5,B5)

A6 = A5.ParentHandle;
B6 = B5.ParentHandle;
IDv(A6,B6)

return
chVector = @(chart)[chart.SpatialSpan, chart.TimeSpan];
chk1 = zeros(atlasB.Size, 4);
for j = 1:atlasB.Size
    chk1(j,:) = chVector(atlasB.Chart(j));
end

for j = 1:atlasB.Size-1
    chV1 = chk1(j,:); % first chart vector
    V1 = sum(abs(chV1 - chk1(j+1:end,:)), 2); % compare with rest of the charts
    if any(V1 < 1e-10)
        Ij = j + find(V1 < 1e-10);
        fprintf('Duplicate Chart: %d and %d \n', [j,Ij])
    end
end




