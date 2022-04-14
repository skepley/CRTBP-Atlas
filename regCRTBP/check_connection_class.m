%CHECK_CONNECTION_CLASS - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
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
% load('test_connection_class.mat')

close all
figure
hold on

% plot primaries 
primarySize = 500;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','b')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')




connections = [sols{:}];

for jConn = 1:size(connections,2)
    cData = connections(:,jConn);
    
    chart1 = cData{1}; 
    chart2 = cData{2};
    localIntersection = cat(1, cData{3}(1:2).',  [cData{3}(3), 0]);
    globalIntersection = cat(1, chart1.local2global(localIntersection(1,:)), chart2.local2global(localIntersection(2,:)))
    phaseIntersection = cell2mat(chart1.eval(localIntersection(1,:))); % A point in R^n which lies on the connection
    
    % swap to F0 field
    switch chart1.RegType
        case 0 % do nothing
        case 1
            phaseIntersection = CRTBP2reg(phaseIntersection(1:4), chart1.Parameter(1), -1);
        case 2
            phaseIntersection = CRTBP2reg(phaseIntersection(1:4), chart1.Parameter(1), -2);
    end % switch
    
    phasePoint = phaseIntersection(1:6);
    flightTime = abs(diff(globalIntersection(:,2)));
    tSpan = linspace(min(globalIntersection(:,2)), max(globalIntersection(:,2)), 100);
    % ob1 = mid(atlas1.orbit(globalIntersection(1,1), tSpan))
    % ob2 = mid(atlas2.orbit(globalIntersection(2,1), tSpan))
    
    
    %  plot connections
    F0 = @(t,x)rk45regvectorfield(t,x,mu,0);
    plot_orbit(@(t,x)-F0(t,x), [0, -.99*globalIntersection(1,2)], phasePoint, [1,3,2], 'PlotOptions', {'Color', [1,0,0]})
    plot_orbit(F0, [0, .99*globalIntersection(2,2)], phasePoint, [1,3,2], 'PlotOptions', {'Color', [0,1,0]})
    scatter3(phasePoint(1), phasePoint(3), phasePoint(2), 'b*')
   
end

% plot primaries 
primarySize = 500;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','b')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')

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
chk1 = zeros(atlas1.Size, 4);
for j = 1:atlas1.Size
    chk1(j,:) = chVector(atlas1.Chart(j));
end

for j = 1:atlas1.Size-1
    chV1 = chk1(j,:); % first chart vector
    V1 = sum(abs(chV1 - chk1(j+1:end,:)), 2); % compare with rest of the charts
    if any(V1 < 1e-10)
        Ij = j + find(V1 < 1e-10);
        fprintf('Duplicate Chart: %d and %d \n', [j,Ij])
    end
end




