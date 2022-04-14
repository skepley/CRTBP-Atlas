%ATLAS_REGCRTBP - testing for RegCRTBPChart atlas construction
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME
 
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 19-Mar-2019; Last revision: 19-Mar-2019

% ================================================== SET UP SYSTEM PATH PROPERTIES ==================================================
whichPC = 'mac';
switch whichPC
    case 'mac' % macbook
        computerPath = '/Users/sk2011/';
    case 'lenovo' % ubuntu
        computerPath = '/home/shane/';
end
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
clear all
clc

% ================================================== SET PARAMETERS FOR THIS ATLAS ==================================================
basis = 'Taylor';
initialTime = 0;
mu = .5;
C = 3.10575;
parameter = [mu;C];
isValid = false;

%% set up some initial data and advect manifold
thetaMid = pi/2; % midpoint of sector
thetaRadius = pi; % radius of sector. setting this to pi gives the entire circle
truncation = [50,15]; 
N = truncation(2);
tau = .1; % initial time re-scaling

%% set up initial data for boundary charts
xInitial = zeros(1,N);
yInitial = zeros(1,N);
rInitial = [1, zeros(1,N-1)];
initCircleRadius = sqrt(8*mu);
expCoeff = exp(1i*thetaMid)*initCircleRadius*(1i*thetaRadius).^(0:N-1)./factorial(0:N-1); % power series for exp(i*pi*s) for s in [-1,1] parameterizes a circle
pInitial = real(expCoeff);
qInitial = imag(expCoeff);
initialData = [xInitial;pInitial;yInitial;qInitial;rInitial]; % Taylor coefficients for initial data

% change to interval arithmetic for testing rigorous solver
if isValid
    initialData = intval(initialData);
end

if exist('bdbackward')
    bdbackward(end+1) = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, 'InitialScaling', -tau, 'boundary', true);
else
    bdbackward = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, 'InitialScaling', -tau, 'boundary', true);
end

if exist('bdforward')
    bdforward(end+1) = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, 'InitialScaling', tau, 'boundary', true);
else
    bdforward = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, 'InitialScaling', tau, 'boundary', true);
end

%% compute atlases for collision manifolds
forwardTau = .33;
forwardAtlas = Atlas(bdforward, tau, forwardTau, @boundarycheck, @advectioncheck);
while ~isempty(forwardAtlas.LeafStack)
    forwardAtlas.growboundary()
end

backwardTau = -.33;
backwardAtlas = Atlas(bdbackward, -tau, backwardTau, @boundarycheck, @advectioncheck);
while ~isempty(backwardAtlas.LeafStack)
    backwardAtlas.growboundary()
end



%% ================================================== PLOT ==================================================
close all
figure;
hold on
tValue = linspace(0,tau,50)';
tRelative = linspace(0, 1, 25)';
sRelative = linspace(-1, 1, 25)';
[S,T] = meshgrid(sRelative, tRelative);
evalData = [reshape(S, [],1), reshape(T, [], 1)];

% plot forward 
% interior and initial boundary
for j = 1:forwardAtlas.Size
    jChart = forwardAtlas.Chart(j);
    [x,p,y,q,r] = jChart.eval(evalData); % Taylor integrator
    if jChart.Generation == 0
        [x0,p0,y0,q0,z0] = jChart.eval([sRelative, zeros(size(sRelative))]);
        plot3(mid(x0),mid(y0),mid(p0), 'k', 'LineWidth',3)
    end
    scatter3(mid(x),mid(y), mid(p), 5,'filled', 'b')
end

% plot terminal boundaries
termChart = RegCRTBPChart;
for jChart = 1:length(forwardAtlas.CrashStack)
    parentChart = forwardAtlas.CrashStack(jChart).ParentHandle;
    termChart(jChart) = parentChart;
    [x1,p1,y1,q1,z1] = parentChart.eval([sRelative, ones(size(sRelative))]);
    plot3(mid(x1),mid(y1),mid(p1), 'r', 'LineWidth',3)
end

% plot bounding circle
bdtheta = linspace(0,2*pi,500);
xb = 0.5*cos(bdtheta);
yb = 0.5*sin(bdtheta);
pb = zeros(size(bdtheta));
plot3(xb,yb,pb,'y','LineWidth',3)

% % plot backward
% % interior and initial boundary
% for j = 1:backwardAtlas.Size
%     jChart = backwardAtlas.Chart(j);
%     [x,p,y,q,r] = jChart.eval(evalData); % Taylor integrator
%     if jChart.Generation == 0
%         [x0,p0,y0,q0,z0] = jChart.eval([sRelative, zeros(size(sRelative))]);
%         plot3(mid(x0),mid(y0),mid(p0), 'k', 'LineWidth',3)
%     end
%     scatter3(mid(x),mid(y), mid(p), 5,'filled','g')
% end


% % plot terminal boundaries
% termChart = RegCRTBPChart;
% for jChart = 1:length(backwardAtlas.CrashStack)
%     parentChart = backwardAtlas.CrashStack(jChart).ParentHandle;
%     termChart(jChart) = parentChart;
%     [x1,p1,y1,q1,z1] = parentChart.eval([sRelative, ones(size(sRelative))]);
%     plot3(mid(x1),mid(y1),mid(p1), 'r', 'LineWidth',3)  
% end


% plot orbits
tf = linspace(0,forwardTau,100);
for s = [-1,1]
    ob = mid(forwardAtlas.orbit(s, tf));
    plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',3)
end

% tb = linspace(0,backwardTau,100);
% ob = mid(backwardAtlas.orbit(0, tb));
% plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',3)

view(194,12)
dealfig()
jobsdone
return

save('collision_manifolds')

%% MINE FOR INTERSECTIONS
sols = {};
fGenIdx = [forwardAtlas.Chart.Generation];
bGenIdx = [backwardAtlas.Chart.Generation];
fGen = 0; % forward generation index
bGen = 0;  % backward generation index
jbGen = backwardAtlas.Chart(bGenIdx == bGen);
nb = length(jbGen);

while fGen <= max(fGenIdx) 
    fGen = fGen + 1;
    jfGen = forwardAtlas.Chart(fGenIdx == fGen);
    disp([fGen,bGen])
    nf = length(jfGen);
    
    % loop over current forward gen
    for ii = 1:nf
        for jj = 1:nb
            jkSols = check4intersection(jfGen(ii), jbGen(jj), 5);
            if ~isempty(jkSols)
                sols{end+1} = {jfGen(ii); jbGen(jj); jkSols};
            end
        end
    end
    
    bGen = bGen + 1;
    jbGen = backwardAtlas.Chart(bGenIdx == bGen);
    disp([fGen,bGen])
    nb = length(jbGen);
    
    % loop over backward gen
    for ii = 1:nf
        for jj = 1:nb
            jkSols = check4intersection(jfGen(ii), jbGen(jj), 5);
            if ~isempty(jkSols)
                sols{end+1} = {jfGen(ii); jbGen(jj); jkSols};
            end
        end
    end
    
end

% save('self_collisions')

return
%% plot self collisions
% load('self_collisions')
clc
% ob = forwardAtlas.Chart(end);
% while ~isempty(ob.ParentHandle)
%     dd = ob.Coordinate.decay;
%     abs(dd(end,:))
%     ob = ob.ParentHandle;
% end


clc
close all
figure
hold on
% plot initial data
for chart = forwardAtlas.Chart
    if isempty(chart.ParentHandle)
        [x0,p0,y0,q0,z0] = chart.eval([sRelative, zeros(size(sRelative))]);
        plot3(mid(x0),mid(y0),mid(p0), 'k', 'LineWidth',3)
    end
end


tForward = linspace(0,.7,500);
tBackward = linspace(0,-.7,500);

% c1 = sols{1}{1}; % first collision chart
% c2 = sols{1}{2};
% [s1,~] = c1.local2global(sols{1}{3}(1),1)
% [s2,~] = c2.local2global(sols{1}{3}(3),1)
% [x,p,y,q,r] = c1.eval(evalData); 
% scatter3(mid(x),mid(y), mid(p), 5,'filled', 'b')
% [x,p,y,q,r] = c2.eval(evalData); 
% scatter3(mid(x),mid(y), mid(p), 5,'filled', 'g')


obf = mid(forwardAtlas.orbit(s1,tForward));
% obb = mid(backwardAtlas.orbit(s2,tBackward));
plot3(obf(:,1),obf(:,3),obf(:,2),'b','LineWidth',1)
% plot3(obb(:,1),obb(:,3),obb(:,2),'g','LineWidth',1)
view(21,70)

% u = obf(:,1);
% p = obf(:,2);
% v = obf(:,3); 
% q = obf(:,4);
% energy = H(u,p,v,q);
% figure
% plot(energy)
%%
tForward = linspace(0,.4,500);
k1 = @(u,v)u.^2 + v.^2;
k2 = @(u,v)u.^2 - v.^2;
U = @(u,v)2*k1(u,v).^3 + 4*mu*k1(u,v).*k2(u,v) + 2*(mu-C).*k1(u,v) + 4*(1-mu) + 4*mu*k1(u,v)./sqrt(k1(u,v).^2 + 1 + 2*k2(u,v));
H = @(u,p,v,q)p.^2 + q.^2 - 2*U(u,v);
close all
figure
hold on
for s = linspace(-1,1,20)
    obf = mid(forwardAtlas.orbit(s,tForward));
    u = obf(:,1);
    p = obf(:,2);
    v = obf(:,3);
    q = obf(:,4);
    r = obf(:,5);
    plot(tForward(1:length(u)), regenergy([u,p,v,q,r],mu,C));
end
legend()
% set(gca,'YLim',[-2,2])



