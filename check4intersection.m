function [isTrue, sols] = check4intersection(chart1, chart2, nInitialCondition, varargin)
%CHECK4INTERSECTION - Use Newton's method to check for an intersection between two regCRTBPCharts
%
%   Syntax:
%       sols = CHECK4INTERSECTION(chart1, chart2, numIC)
%
%   Inputs:
%       chart1 - backward time chart from a stable manifold of the form P(s,t)
%       chart2 - forward time chart from an unstable manifold of the form Q(s,t)
%       nInitialCondition - Description
%
%   Outputs:
%       sols - a vector of the form x = (s1, t1, s2) such that Q(s2, t0) = P(s1, t1) in the 3 coordinates specified.

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 21-Mar-2019; Last revision: 10-Mar-2021

%% parse input
p = inputParser;
addRequired(p,'chart1')
addRequired(p,'chart2')
addRequired(p,'numIC')
addParameter(p,'dropCoordinate', 4)
addParameter(p,'t0', 0);
addParameter(p,'NumIterate', 5);

parse(p,chart1, chart2, nInitialCondition, varargin{:})
dropCoordinate = p.Results.dropCoordinate;
t0 = p.Results.t0;
numIterate = p.Results.NumIterate; % look for connections of stable manifold up to tau = timeInterval(1) and unstable manifold up to tau = timeInterval(2)


isTrue = false;  % return false by default even when no solution is found
% First try to rule out intersection using ell-1 box method. If this fails it tries to find an intersection via Newton iteration.
% if isequal(chart1.RegType, chart2.RegType)  % this is only implemented so far for same regularization types

if ~intersectellonebox(chart1, chart2) % intersectellonebox returns false if intersection has been ruled out
    sols = [];
    return
    % elseif ~isequal(chart1.RegType, chart2.RegType) % charts are not in same coordinates
    %
    %     if isequal(chart1.RegType, 0) || isequal(chart1.RegType, 0) % exactly 1 chart is in F0 coordinates
    %         if isequal(chart1.RegType, 0) % chart1 is in F0 coordinates
    %             chart1 = chart1.hotswap(chart2.RegType); % change chart 1 from F0 into F1 or F2 coordinates
    %         else % chart2 is in F0 coordinates
    %             chart2 = chart2.hotswap(chart1.RegType); % change chart 2 from F0 into F1 or F2 coordinates
    %         end
    %
    %     else % charts are in F1 and F2 coordinates
    %         chart1 = chart1.hotswap(0); % change chart 1 to F0 coordinates
    %         chart2 = chart2.hotswap(0); % change chart 2 to F0 coordinates
    %     end
end
% end

%% ell-1 box method did not rule out an intersection so move on to Newton iteration
% set projection map
coordIdx = 1:4;
projmap = @(u)u(coordIdx ~= dropCoordinate, :);  % projection onto only 3 rows. Apply to column vectors or matrices.

%% Fix the time variable and define Q (unstable chart) and its partial derivatves as functions of a single variable, r.
Qr = chart2.Coordinate.fixtime(t0); % collapse onto t = t0 and evaluate
Qx = polynom(flip(Qr(1).Coefficient),'r');
Qp = polynom(flip(Qr(2).Coefficient),'r');
Qy = polynom(flip(Qr(3).Coefficient),'r');
Qq = polynom(flip(Qr(4).Coefficient),'r');
Q = @(r)mid([polyval(Qx,r); polyval(Qp,r); polyval(Qy,r); polyval(Qq, r)]);

dQx = pderiv(Qx);
dQp = pderiv(Qp);
dQy = pderiv(Qy);
dQq = pderiv(Qq);
dQdr = @(r)mid([polyval(dQx,r); polyval(dQp,r); polyval(dQy,r); polyval(dQq, r)]);

%% Define P (stable chart) and its partial derivatives in each direction evaluated as functions of two variables, (s,t).
Pst = chart1.Coordinate;
Px = intlabpoly(Pst(1));
Pp = intlabpoly(Pst(2));
Py = intlabpoly(Pst(3));
Pq = intlabpoly(Pst(4));
P = @(s,t)mid([polyval(Px,[s,t]); polyval(Pp,[s,t]); polyval(Py,[s,t]); polyval(Pq, [s,t])]);

% dPdt
dtPx = pderiv(Px,'t');
dtPp = pderiv(Pp,'t');
dtPy = pderiv(Py,'t');
dtPq = pderiv(Pq, 't');
dPdt = @(s,t)mid([polyval(dtPx,[s,t]); polyval(dtPp,[s,t]); polyval(dtPy,[s,t]); polyval(dtPq, [s,t])]);

% dPds
dsPx = pderiv(Px,'s');
dsPp = pderiv(Pp,'s');
dsPy = pderiv(Py,'s');
dsPq = pderiv(Pq, 's');
dPds = @(s,t)mid([polyval(dsPx,[s,t]); polyval(dsPp,[s,t]); polyval(dsPy,[s,t]); polyval(dsPq, [s,t])]);

%% Set up Zero finding problem and Newton map and iterate

% define zero finding map, derivative, and Newton operator
PQregType = [chart1.RegType, chart2.RegType];
mu = chart1.Parameter(1);
[F, DF, N] = setnewtonmap(mu, P, Q, dPds, dPdt, dQdr, projmap, PQregType);
warning('off','MATLAB:nearlySingularMatrix') % turn off warning for ill conditioned jacobian
sols = []; % initialize solution vector
iter = 1;
while iter <= nInitialCondition && isempty(sols) % no solution found yet
    x0 = 2*rand(3,1) -1; % initial condition in [-1,1]^3
    jj = 1;
    while jj <= numIterate && norm(x0,inf) <= 10
        jj = jj+1;
        x0 = N(x0); % (s1,t1,s2) such that P(s1,t1) (chart 1) intersects Q(s2,t0) (chart2)
    end
    y = F(x0);
    
    % Verify the following: 1. F(x) < 1e-5, 2. DF(x0) is non-singular, 3. x0 in [-1,1]x[0,1]x[-1,1]
    if norm(y) < 1e-5 && norm(x0,inf) <= 1 && x0(2) >= 0 && rank(DF(x0)) > 2 %
        for jj = 1:10 % 10 more rounds of Newton to improve approximation
            x0 = N(x0);
        end
        sols = x0;
        
        % evaluate both charts at candidate intersection point and check if it is a true solution or pseudo-solution
        c1 = P(x0(1), x0(2));
        c2 = Q(x0(3));
        dTrue = c1 - c2; % exactly zero at a true intersection
        dPseudo = c1 + c2; % exactly zero at a pseudo-intersection
        isTrue = abs(dTrue(4)) < abs(dPseudo(4)); % a candidate solution if it is closer satisfying the true condition
        if isTrue
            disp('true-solution')
        else
            disp('pseudo-solution')
        end
    end
    iter = iter + 1; % increment the iteration count
end
end % end check4intersection

function [F, DF, N] = setnewtonmap(mu, P, Q, dPds, dPdt, dQdr, projmap, PQregType)
% Return functions for evaluating a zero finding problem (F), its differential (DF), and the associated Newton iteration map (N)
PQRegTypeCase = 3*PQregType(1) + PQregType(2);  % treat PQregType as a ternary expansion and convert integer case number
trimAD = @(u)u(1:4);  % remove automatic differentiation coordinates from the regularization map implementation
switch PQRegTypeCase
    case 0 % Both charts are in f_0 coordinates so apply the naive zero finding map
        F = @(u)projmap(P(u(1),u(2)) - Q(u(3)));
        DF = @(u)projmap([dPds(u(1),u(2)), dPdt(u(1),u(2)), -dQdr(u(3))]);
        
    case  1 % [0, 1] compose Q with regularization map, u_1 |--> u_0
        regMap = @(u)trimAD(CRTBP2reg(u.', mu, -1)).';
        F = @(u)projmap(P(u(1),u(2)) - regMap(Q(u(3))));
        DF = @(u)projmap(cat(2, [dPds(u(1),u(2)), dPdt(u(1),u(2))],...
            -diffCRTBP2reg(Q(u(3)))*dQdr(u(3))));
        
    case 2  % [0, 2] compose Q with regularization map, u_2 |--> u_0
        regMap = @(u)trimAD(CRTBP2reg(u.', mu, -2)).';
        F = @(u)projmap(P(u(1),u(2)) - regMap(Q(u(3))));
        DF = @(u)projmap(cat(2, [dPds(u(1),u(2)), dPdt(u(1),u(2))],...
            -diffCRTBP2reg(Q(u(3)))*dQdr(u(3))));
        
    case 3  % [1, 0] compose P with regularization map, u_1 |--> u_0
        regMap = @(u)trimAD(CRTBP2reg(u.', mu, -1)).';
        F = @(u)projmap(regMap(P(u(1),u(2))) - Q(u(3)));
        DF = @(u)projmap(cat(2, diffCRTBP2reg(P(u(1), u(2)))*[dPds(u(1),u(2)), dPdt(u(1),u(2))],...
            -dQdr(u(3))));
        
    case 6  % [2, 0]  compose P with regularization map, u_2 |--> u_0
        regMap = @(u)trimAD(CRTBP2reg(u.', mu, -2)).';
        F = @(u)projmap(regMap(P(u(1),u(2))) - Q(u(3)));
        DF = @(u)projmap(cat(2, diffCRTBP2reg(P(u(1), u(2)))*[dPds(u(1),u(2)), dPdt(u(1),u(2))],...
            -dQdr(u(3))));
        
        %     case 5  % [1, 2] compose both P and Q with inverse regularization maps taking them to u_0
        %         PregMap = @(u)trimAD(CRTBP2reg(u.', mu, -1)).';
        %         QregMap = @(u)trimAD(CRTBP2reg(u.', mu, -2)).';
        %         F = @(u)projmap(PregMap(P(u(1),u(2))) - QregMap(Q(u(3))));
        %         DF = @(u)projmap(cat(2, diffCRTBP2reg(P(u(1), u(2)))*[dPds(u(1),u(2)), dPdt(u(1),u(2))],...
        %             -diffCRTBP2reg(Q(u(3)))*dQdr(u(3))));
        
        %     case 7  % [2, 1]  compose both P and Q with inverse regularization maps taking them to u_0
        %         PregMap = @(u)trimAD(CRTBP2reg(u.', mu, -2)).';
        %         QregMap = @(u)trimAD(CRTBP2reg(u.', mu, -1)).';
        %         F = @(u)projmap(PregMap(P(u(1),u(2))) - QregMap(Q(u(3))));
        %         DF = @(u)projmap(cat(2, diffCRTBP2reg(P(u(1), u(2)))*[dPds(u(1),u(2)), dPdt(u(1),u(2))],...
        %             -diffCRTBP2reg(Q(u(3)))*dQdr(u(3))));
        
    case {4, 5, 7, 8} % {[1, 1], [1,2], [2,1], [2, 2]} compose both P and Q with inverse regularization maps taking them to u_0
        PregMap = @(u)trimAD(CRTBP2reg(u.', mu, -PQregType(1))).';
        QregMap = @(u)trimAD(CRTBP2reg(u.', mu, -PQregType(2))).';
        F = @(u)projmap(PregMap(P(u(1),u(2))) - QregMap(Q(u(3))));
        DF = @(u)projmap(cat(2, diffCRTBP2reg(P(u(1), u(2)))*[dPds(u(1),u(2)), dPdt(u(1),u(2))],...
            -diffCRTBP2reg(Q(u(3)))*dQdr(u(3))));
    otherwise
        error('This should not happen')
end
N = @(u)(u - DF(u)\F(u)); % Define the Newton map
end
% Revision History:
%{
10-Mar-2021 - Added implementation for mining intersection in charts of different regularization type.
%}


