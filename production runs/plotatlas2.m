%PLOTATLAS2 - Compute plots of the (un)stable atlases computed in final2.m which must be run first.
%   Plot (un)stable manifolds for the following example:
%       Parameters: equal mass
%       Energy: L4
%       Stable manifold: L4
%       Unstable manifold: L4
%       Time: +/- 7 units
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 06-Feb-2020; Last revision: 08-Mar-2021




%% EST BED FOR THE CRTBPpatch method
close all
figure
hold on

% plot primaries and L4
primarySize = 100;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','b')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')
plot3(L4(1), L4(3), L4(2), 'r*')


plotIdx = [1, 3, 2];
unstableColor = [0, 1, 0];
stableColor = [1, 0, 0];

globalSpace = linspace(-1, 1, 100);
unstableTime = linspace(0, maxTau, 100);
stableTime = linspace(0, minTau, 100);
unstableEvalArray = {globalSpace, unstableTime};
stableEvalArray = {globalSpace, stableTime};



atlasF.patch(unstableEvalArray, 0, 0, plotIdx, 'FaceColor', unstableColor);
atlasB.patch(stableEvalArray, 0, 0, plotIdx, 'FaceColor', stableColor);

return
%%  ================================================== Function block begins here ==================================================
hold on 
for iChart = atlasB.Chart% loop through charts of correct regtype
    if isequal(iChart.RegType, 0)
        iChart.plotboundary([40,40], 0, plotIdx, 'PlotOptions', {'LineWidth', 1, 'Color', [0, 0, 1]})
    end
end

%  plot primaries and L4
primarySize = 500;
scatter3(mu,0,0,(1-mu)*primarySize,'filled','k')
scatter3(mu-1,0,0,mu*primarySize,'filled','k')
plot3(L4(1), L4(3), L4(2), 'r*')
plot(L4(1), L4(3), 'r*')

axis tight
ax = axis;
dealfig()

