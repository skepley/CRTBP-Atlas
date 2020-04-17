function sols = check4intersection(chart1, chart2, nInitialCondition, varargin)
%CHECK4INTERSECTION - Use Newton's method to check for an intersection between two regCRTBPCharts
%
%   CHECK4INTERSECTION() - A more detailed description of the function
%
%   Syntax:
%       sols = CHECK4INTERSECTION(chart1, chart2, numIC)
%
%   Inputs:
%       chart1 - forward time chart from an unstable manifold
%       chart2 - backward time chart from a stable manifold
%       nInitialCondition - Description
%
%   Outputs:
%       sols - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 21-Mar-2019; Last revision: 21-Mar-2019

%% parse input
p = inputParser;
addRequired(p,'chart1')
addRequired(p,'chart2')
addRequired(p,'numIC')
addParameter(p,'KeepCoord',2)
addParameter(p,'t0',0);
addParameter(p,'NumIterate',5);

parse(p,chart1, chart2, nInitialCondition, varargin{:})
keepCoord = p.Results.KeepCoord;
t0 = p.Results.t0;
numIterate = p.Results.NumIterate; % look for connections of stable manifold up to tau = timeInterval(1) and unstable manifold up to tau = timeInterval(2)

% First try to rule out intersection using ell-1 box method. If this fails it tries to find an intersection via Newton iteration.
if ~intersectellonebox(chart1, chart2) % intersectellonebox returns false if intersection has been ruled out
    sols = [];
    return
elseif ~isequal(chart1.RegType, chart2.RegType) % charts are not in same coordinates
    
    if isequal(chart1.RegType, 0) || isequal(chart1.RegType, 0) % exactly 1 chart is in F0 coordinates
        if isequal(chart1.RegType, 0) % chart1 is in F0 coordinates
            chart1 = chart1.hotswap(chart2.RegType); % change chart 1 from F0 into F1 or F2 coordinates
        else % chart2 is in F0 coordinates
            chart2 = chart2.hotswap(chart1.RegType); % change chart 2 from F0 into F1 or F2 coordinates
        end
        
    else % charts are in F1 and F2 coordinates
        chart1 = chart1.hotswap(0); % change chart 1 to F0 coordinates
        chart2 = chart2.hotswap(0); % change chart 2 to F0 coordinates
    end
end

% get dropped velocity coordinate
vCoord = [2,4];
dropCoord = vCoord(vCoord ~= keepCoord);

% fix time for Q and get derivative
Qr = chart2.Coordinate.fixtime(t0); % collapse onto t = t0
Q1 = polynom(flip(Qr(1).Coefficient),'r');
Q2 = polynom(flip(Qr(3).Coefficient),'r');
Q3 = polynom(flip(Qr(keepCoord).Coefficient),'r'); % third coordinate is one of the velocity components
Q = @(r)[polyval(Q1,r);polyval(Q2,r);polyval(Q3,r)];

dQ1 = pderiv(Q1);
dQ2 = pderiv(Q2);
dQ3 = pderiv(Q3);
dQ = @(r)[polyval(dQ1,r);polyval(dQ2,r);polyval(dQ3,r)];

% compute P and partial derivatives of P in each direction
Pst = chart1.Coordinate([1,3,keepCoord]); % throw out one velocity coordinate
P1 = intlabpoly(Pst(1));
P2 = intlabpoly(Pst(2));
P3 = intlabpoly(Pst(3));
P = @(s,t)[polyval(P1,[s,t]);polyval(P2,[s,t]);polyval(P3,[s,t])];

% dPdt
dtP1 = pderiv(P1,'t');
dtP2 = pderiv(P2,'t');
dtP3 = pderiv(P3,'t');
dPdt = @(s,t)[polyval(dtP1,[s,t]);polyval(dtP2,[s,t]);polyval(dtP3,[s,t])];

% dPds
dsP1 = pderiv(P1,'s');
dsP2 = pderiv(P2,'s');
dsP3 = pderiv(P3,'s');
dPds = @(s,t)[polyval(dsP1,[s,t]);polyval(dsP2,[s,t]);polyval(dsP3,[s,t])];

% function to check dropped coordinate on Newton solutions
P4 = intlabpoly(chart1.Coordinate(dropCoord));
Q4 = intlabpoly(chart2.Coordinate(dropCoord));
F4 = @(x)abs(polyval(P4, x(1:2)) - polyval(Q4, x(3:4))); % F4: R^4 ---> R checks if the last coordinate is correct
F4energy = @(x)abs(polyval(P4, x(1:2)) + polyval(Q4, x(3:4))); % make sure false positives are actually in the same energy

% define zero finding map, derivative, and Newton operator
F = @(v)mid(P(v(1),v(2)) - Q(v(3)));
DF = @(v)mid([dPds(v(1),v(2)), dPdt(v(1),v(2)), -dQ(v(3))]);
N = @(x)(x - DF(x)\F(x)); % Newton iteration
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
        for jj = 1:10 % 10 more rounds to Newton to improve approximation
            x0 = N(x0);
        end
        
        % get last coordinate for charts at intersection candidate
        q1 = polyval(P4, x0(1:2)); % P4(s1,t1)
        q2 = polyval(Q4, [x0(3), t0]); % Q4(s2, t0)
        dTrue = q1 - q2; % exactly zero at a true intersection
        dPseudo = q1 + q2; % exactly zero at a pseudo-intersection
        if abs(dTrue) < abs(dPseudo) % a candidate solution if it is closer satisfying the true condition
            sols = x0;
            disp('true-solution')
        else
            disp('pseudo-solution')
        end
    end
    iter = iter + 1; % increment the iteration count
end % end check4intersection

% Revision History:
%{

%}


% ============================ OLD VERSION ============================
%
%         fullSol = zeros(1,4); % initialize full solution with configuration coordinates
%         fullSol(1:3) = x0;
%         fullSol(4) = t0;
%         if F4(fullSol) < 1e-5
%             sols = x0;
%             disp('true-solution')
%         else
%             disp('pseudo-solution')
%
%             if F4energy(fullSol) > 1e-5
%                 disp('something is wrong here')
%             end
%         end


% ============================ OLD OLD VERSION ============================
%         y = F(x0);
%         p0 = chart1.eval(x0(1:2)'); % point of intersection through each chart
%         q0 = chart2.eval([x0(3),t0]);
%         if norm(p0 - q0) < 1e-3 % && 1e-2 < abs(p0(4) + q0(4)) % intersections are close, and missing velocity has correct sign and is not nearly zero.
%             if norm(y) < 1e-10 && norm(x0,inf) <= 1 && x0(2) <= 1e-12 % final check to make sure output is a zero in [-1,1]x[-1,0]x[-1,1]
%                 sols = x0;
%             else
%                 disp('Final connection check failed. Something is wrong')
%                 sols = [];
%             end
%         elseif norm(p0 - q0) < 1e-3 && abs(p0(4)) < 1e-3 && isequal(keepCoord,2) % intersection is close, v1 was kept but is nearly zero at intersection
%             fprintf('sign of v1 is indeterminate, switching to v2 \n')
%             sols = checkconnection(chart1,chart2,nInitialCondition,'keepCoord',4,'t0',t0,'NumIterate',numIterate);
%         elseif norm(p0 - q0) < 1e-3 && abs(p0(4)) < 1e-3 && isequal(keepCoord,4) % intersection is close, both velocities nearly zero at intersection
%             fprintf('sign of v1 and v2 are indeterminate. These charts must be followed through the atlas \n')
%             sols = [];

