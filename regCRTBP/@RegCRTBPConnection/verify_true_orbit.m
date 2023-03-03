function scores = verify_true_orbit(obj)
%VERIFY_TRUE_ORBIT - Return some heuristics for guessing whether or not a connection is true or pseudo.
%
%   Syntax:
%       [tangentScore, rk45_Score] = VERIFY_TRUE_ORBIT(obj)
%
%   Inputs:
%       obj - RegCRTBPConnection
%
%   Outputs:
%       tangentScore - cos(theta) where theta is the angle between the vector field evaluated at the intersection point using the 
%           stable and unstable charts
%       rk45_score - The sum of distances between the time-tau map for both initial conditions with tau = +/- 0.5 units
%
%   Subfunctions: none
%   Classes required: RegCRTBPConnection
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 24-Jan-2023;


score1 = tangentScore(obj);
score2 = rk45_Score(obj);
scores = [score1, score2];

end % end verify_true_orbit

function score = tangentScore(C)
% score connection based on how close to tangent the vector field is at the two intersection images
mu = C.Parameter;
if isequal(C.StableChart.RegType, C.UnstableChart.RegType)
    regType = C.StableChart.RegType;
else
    error('This is not implemented yet')
end

F = @(t,x)rk45regvectorfield(t, x, mu, regType);  % rk45 field to test against
stableInitialData = cell2mat(C.StableChart.eval(C.LocalIntersection(1:2).'));
unstableInitialData = cell2mat(C.UnstableChart.eval(C.LocalIntersection(3:4).'));
Fs = F(0, stableInitialData);
Fu = F(0, unstableInitialData);
score = dot(Fs, Fu)/(norm(Fs)*norm(Fu));  % return cos(theta) where theta is the angle between Fu and Fs
end

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