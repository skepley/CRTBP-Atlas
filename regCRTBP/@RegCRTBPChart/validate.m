function varargout = validate(obj, varargin)
%VALIDATE - Computes validation of a single chart for Lorenz integrator (with bootstrapping)

%   Description: The vector field is given by
%           VECTOR FIELD:
%           xdot = p
%           pdot = 8(x^2 + y^2)q + 12x(x^2 + y^2)^2 + 16mu*x^3 + 4(mu - C)x + 8mu*(x^2 - 3y^2 + 1)xr^3
%           ydot = q
%           qdot = ?8(x^2 + y^2)p + 12y(x^2 + y^2)^2 - 16mu*y^3 + 4(mu - C)y + 8mu*(-y^2 + 3x^2 + 1)yr^3
%           rdot = -2r^3((x^2 + y^2)(xp + yq) + (xp - yq))
%
%   Inputs:
%       obj - instance of the RegCRTBPChart class
%       tau - size of timestep (also equal to the vector field rescaling)
%
%   Subfunctions: none
%   Classes required: @Chart, @Scalar
%   Other m-files required: solveradiipoly
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 07-Mar-2019; Last revision: 18-Mar-2019

if ~strcmp(obj.NumericalClass, 'intval')
    warning('validation is not rigorous without interval coefficient Scalars')
end

if nargin > 1
    printYZBounds = varargin{1}; % print Y, Z and error bounds to screen
else
    printYZBounds = false;
end

%% set up intval parameters, scaling, and truncation
mu = obj.Parameter(1); % these are the same numerical class as the chart
C = obj.Parameter(2);
M = obj.Truncation(1); % number of temporal terms
tau = intval(obj.Tau); % vector field scaling

% define basis
basis = obj.Coordinate(1).Basis;
if ~strcmp(basis, 'Taylor')
    error('validate - not implemented for non-Taylor basis')
end

%% YO BOUNDS (DEFECT)
Yscaling = tau;
fulltimes = @(a,b)mtimes(a,b, 'Full');

% generate base elements of order 1 for approximate solution
x = obj.Coordinate(1);
p = obj.Coordinate(2);
y = obj.Coordinate(3);
q = obj.Coordinate(4);
r = obj.Coordinate(5);

% generate base monomials of order 2 (no truncation)
xx = fulltimes(obj.Coordinate(1), obj.Coordinate(1)); % x^2
yy = fulltimes(obj.Coordinate(3), obj.Coordinate(3)); % y^2
xx_plus_yy = xx + yy; % x^2 + y^2
rr = fulltimes(obj.Coordinate(5), obj.Coordinate(5)); % r^2
xp = fulltimes(obj.Coordinate(1), obj.Coordinate(2)); % x*p
yq = fulltimes(obj.Coordinate(3), obj.Coordinate(4)); % y*q
xp_plus_yq = xp + yq; % x*p + y*q
xp_minus_yq = xp - yq; % x*p - y*q
h21 = xx - 3*yy + 1; % x^2 - 3y^2 + 1
h22 = 3*xx - yy + 1; % -y^2 + 3x^2 + 1

% generate base monomials of order 3 (no truncation)
xxx = fulltimes(obj.Coordinate(1), xx); % x^3
yyy = fulltimes(obj.Coordinate(3), yy); % y^3
rrr = fulltimes(obj.Coordinate(5), rr); % r^3
h31 = fulltimes(xx_plus_yy, obj.Coordinate(4)); % (x^2 + y^2)q
h32 = fulltimes(xx_plus_yy, obj.Coordinate(2)); % (x^2 + y^2)p

% generate base monomials of order 4 (no truncation)
xrrr = fulltimes(obj.Coordinate(1), rrr); % xr^3
yrrr = fulltimes(obj.Coordinate(3), rrr); % yr^3
h41 = fulltimes(xx_plus_yy, xx_plus_yy); % (x^2 + y^2)^2
h42 = fulltimes(xx_plus_yy, xp_plus_yq); % (x^2 + y^2)(xp + pq)
h43 = h42 - xp_minus_yq; % (x^2 + y^2)(xp + pq) - (xp - pq)

% generate base monomials of order 5 (no truncation)
h51 = fulltimes(h41, obj.Coordinate(1)); % x(x^2 + y^2)^2
h52 = fulltimes(h41, obj.Coordinate(3)); % y(x^2 + y^2)^2

% generate base monomials of order 6 (no truncation)
h61 = fulltimes(h21, xrrr); % (x^2 - 3y^2 + 1)xr^3
h62 = fulltimes(h22, yrrr); % (-y^2 + 3x^2 + 1)yr^3

% generate base monomials of order 7 (no truncation)
h71 = fulltimes(rrr, h43); %  r^3((x^2 + y^2)(xp + pq) - (xp - pq))


% define tail projection operator
if isequal(obj.Dimension(1), 1) % initial data is a point.
    tailProjection = @(scalar)Scalar(scalar.Coefficient(M+1:end), basis); % pi^Inf
elseif isequal(obj.Dimension(1), 2) % initial data is an arc
    N = obj.Truncation(2);
    tailProjection = @(scalar)Scalar(scalar.Coefficient(M+1:end, N+1:end), basis); % pi^Inf
else
    error('validate - Y0 bounds for higher dimensional initial data not implemented')
end

% rigorous bound for on norm for T(0,0,0)
Y1 = norm(tailProjection(int(p)));
Y2 = norm(tailProjection(int(8*h31 + 12*h51 + xxx*(16*mu) + x*(4*(mu-C)) + h61*(8*mu))));
Y3 = norm(tailProjection(int(q)));
Y4 = norm(tailProjection(int(-8*h32 + 12*h52 - yyy*(16*mu) + y*(4*(mu-C)) + h62*(8*mu))));
Y5 = norm(tailProjection(int(-2*h71)));
Y0 = Yscaling*max([Y1,Y2,Y3,Y4,Y5]) + obj.ErrorBound;

%% Z1 bounds (DERIVATIVE)
Zscaling = Yscaling/(M+1);
normVector = obj.Coordinate.norm(); % (||x||, ||p||, ||y||, ||q||, ||r||) in R^5
nX = normVector(1);
nP = normVector(2);
nY = normVector(3);
nQ = normVector(4);
nR = normVector(5);

n1 = xx_plus_yy.norm(); % ||x^2 + y^2||
n2 = h21.norm(); % ||x^2 - 3y^2 + 1||
n3 = h22.norm(); % ||y^2 - 3x^2 + 1||
n4 = xp_plus_yq.norm(); % ||xp + yq||
n5 = xp_minus_yq.norm(); % ||xp - yq||

DT1 = intval([0,1,0,0,0]);
DT2 = [16*nQ*nX + 12*(n1^2 + 4*nX^2*n1) + 48*mu*nX^2 + 4*(mu-C) + 8*mu*nR^3*(n2 + 2*nX^2),... % DT11
    intval(0), 16*nQ*nY + 48*nX*nY*n1 + 48*mu*nR^3*nX*nY,... % DT12, DT13
    8*n1, 24*mu*nR^2*nX*n2]; % DT14, DT15
DT3 = intval([0,0,0,1,0]);
DT4 = [16*nP*nX + 48*nX*nY*n1 + 48*mu*nR^3*nX*nY, 8*n1,... % DT31, DT32
    16*nP*nY + 12*(n1^2 + 4*nY^2*n1) + 48*mu*nY^2 + 4*(mu-C) + 8*mu*nR^3*(n3 + 2*nY^2),... % DT33
    intval(0), 24*mu*nR^2*nY*n3]; % DT34, DT35
DT5 = cat(2,2*nR^3*[2*nX*n4 + nP*(n1 + 1), nX*(n1 + 1), 2*nY*n4 + nQ*(n1 + 1), nY*(n1 + 1)], 6*nR^2*(n1*n4 + n5)); % DT5
ZMatrix = [DT1; DT2; DT3; DT4; DT5];
Z1 = Zscaling*norm(ZMatrix, Inf);

%% Z2 bounds (DERIVATIVE)
rStar = intval(10)*Y0;
nXStar = nX + rStar;
nPStar = nP + rStar;
nYStar = nY + rStar;
nQStar = nQ + rStar;
nRStar = nR + rStar;
n1Star = nXStar^2 + nYStar^2; % (||x||+r)^2 + (||y||+r)^2
n2Star = n1Star + 2*nYStar^2 + intval(1); % (||x||+r)^2 + 3(||y||+r)^2 + 1
n3Star = n1Star + 2*nXStar^2 + intval(1); % (||y||+r)^2 + 3(||x||+r)^2 + 1
n4Star = nXStar*nPStar + nYStar*nQStar; % (||x||+r)(||p||+r) + (||y||+r)(||q||+r) 
n5Star = nXStar*nQStar + nYStar*nPStar; % (||x||+r)(||q||+r) + (||y||+r)(||p||+r) 


H1 = intval(zeros(5,5)); % Hessian of T_1 evaluated at norm (||a|| + rStar)
H2 = intval([16*nQStar + 24*(n1Star + 4*(2*nXStar^3 + nXStar*nYStar^2) + 4*mu*nXStar + 2*mu*nRStar^3*nXStar),... % h11
    0, 24*(2*n1Star*nYStar + 4*nXStar^2*nYStar + 4*mu*nRStar^3*nYStar), 16*nXStar, 24*mu*nRStar^2*(3*n1Star + 1);... % h12,h13,h14,h15
    0, 0, 0, 0, 0;... % row 2
    48*nYStar*(3*nXStar^2 + nYStar^2 + mu*nRStar^3), intval(0), 16*nQStar + 48*nXStar*(nXStar^2 + 3*nYStar^2 + mu*nR^3),... % h31,h32,h33
    16*nYStar, 144*mu*nXStar*nYStar*nRStar^2;... % h34, h35
    16*nXStar, 0, 16*nYStar, 0, 0;... % row 4
    24*mu*(3*n1Star + 1)*nRStar^2, 0, 144*mu*nXStar*nYStar, 0, 48*mu*nXStar*(nXStar^2 + 3*nYStar + 1)*nRStar]); % row 5
H3 = intval(zeros(5,5));
H4 = intval([16*nPStar + 24*(n1Star + 4*(2*nYStar^3 + nXStar^2*nYStar) + 4*mu*nYStar + 2*mu*nRStar^3*nYStar),... % h11
    0, 24*(2*n1Star*nXStar + 4*nXStar*nYStar^2 + 4*mu*nRStar^3*nXStar), 16*nYStar, 24*mu*nRStar^2*(3*n1Star + 1);... % h12,h13,h14,h15
    0, 0, 0, 0, 0;... % row 2
    48*nXStar*(nXStar^2 + 3*nYStar^2 + mu*nRStar^3), intval(0), 16*nPStar + 48*nYStar*(3*nXStar^2 + nYStar^2 + mu*nR^3),... % h31,h32,h33
    16*nXStar, 144*mu*nXStar*nYStar*nRStar^2;... % h34, h35
    16*nYStar, 0, 16*nXStar, 0, 0;... % row 4
    24*mu*(3*n1Star + 1)*nRStar^2, 0, 144*mu*nXStar*nYStar, 0, 48*mu*nYStar*(3*nXStar^2 + nYStar + 1)*nRStar]); % row 5
H5 = intval([4*nRStar^3*(3*nXStar*nPStar + nYStar*nQStar), 2*nRStar^3*(3*nXStar^2 + nYStar^2 + 1), 4*nRStar^3*n5Star,... % h11, h12, h13
    4*nRStar^3*nXStar*nYStar, 6*nRStar^2*(2*nXStar*n4Star + nPStar*(n1Star + 1));... % h14, h15
    2*nRStar^3*(3*nXStar^2 + nYStar^2 + 1), 0, 4*nRStar^3*nXStar*nYStar, 0, 6*nRStar^2*nXStar*(n1Star + 1);...% row 2
    4*nRStar^3*n5Star, 4*nRStar^3*nXStar*nYStar, 4*nRStar^3*n4Star, 2*nRStar^3*(nXStar^2 + 3*nYStar^2 + 1), 6*nRStar^2*(2*nYStar*n4Star + nQStar*(n1Star + 1));... % row 3
    4*nRStar^3*nXStar*nYStar, 0, 2*nRStar^3*(nXStar^2 + 3*nYStar^2), 0, 6*nRStar^2*nYStar*(n1Star + 1);... % row 4
    6*nRStar^2*(2*nXStar*n4Star + nPStar*(n1Star + 1)), 6*nRStar^2*nXStar*(n1Star + 1), 6*nRStar^2*(2*nYStar*n1Star + nQStar*(n1Star + 1)),... %h51,h52,h53
    6*nRStar^2*nYStar*(n1Star + 1), 12*nRStar*(n1Star*(n4Star + 1))]);


% Each element of ZVector is the norm of T(aBar + u) evaluated with the real scalar bound ||u|| = rStar. 
ZVector(1) = norm(H1, Inf); % ||D^2_T1(x,y)|| 
ZVector(2) = norm(H2, Inf); % ||D^2_T1(x,y)||
ZVector(3) = norm(H3, Inf); % ||D^2_T3(x,y)||
ZVector(4) = norm(H4, Inf); % ||D^2_T4(x,y)||
ZVector(5) = norm(H5, Inf); % ||D^2_T5(x,y)||
Z2 = Zscaling*norm(ZVector, Inf);

%% Solve radii polynomial
%  Find interval of negativity for radii polynomial: p(r) = Z2r^2 - (1-Z1)r + Y0

if any([Y0, Z1, Z2] == Inf)
    warning('Y0, Z1, or Z2 is infinite')
    maxError = Inf;
else
    r1 = solveradiipoly([Z2,-(1-Z1),Y0], printYZBounds);
    maxError = sup(r1);
end

rStarFlag = false;
if exist('rStar') && rStar < maxError
    rStarFlag = true;
    warning('r* < r')
end

if ~isreal(maxError) || maxError < 0 || rStarFlag
    obj.ErrorBound = Inf; % validation failed, set error to Inf
    obj.Crash = 1; % stop further integration
    warning('validation failed')
else
    obj.ErrorBound = obj.InitialError + maxError; % error from this step + previous propagation
end

if printYZBounds
    fprintf('r: %.3g \n', obj.ErrorBound)
    fprintf('Y0: %.3g \n', mid(Y0))
    fprintf('Z1: %.3g \n', mid(Z1))
    fprintf('Z2: %.3g \n', mid(Z2))
end

switch nargout
    case 0
        return
    case 1
        varargout{1} = obj.ErrorBound;
    case 2
        varargout{1} = obj.ErrorBound;
        varargout{2} = [Y0, Z1, Z2];
end % end switch
end % end validate
