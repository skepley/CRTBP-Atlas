function varargout = validate(obj, varargin)
%VALIDATE - Computes validation of a single chart for Lorenz integrator (with bootstrapping)

%   Description: The vector field is given by
%           VECTOR FIELD:
%           xdot = p
%           pdot = 2q + mu1(x-mu2)(1-r1^3) + mu2(x+mu1)(1-r2^3)
%           ydot = q
%           qdot = -2p + mu1y(1-r1^3) + mu2y(1-r2^3)
%           r1dot = -r1^3((x-mu2)p + yq)
%           r2dot = -r2^3((x-mu2)p + yq)
%
%   Inputs:
%       obj - instance of the CRTBPChart class
%       tau - size of timestep (also equal to the vector field rescaling)
%
%   Subfunctions: none
%   Classes required: @Chart, @Scalar, @CRTBPChart
%   Other m-files required: solveradiipoly
%   MAT-files required: none
%
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 20-Mar-2019; Last revision: 20-Mar-2019

if ~strcmp(obj.NumericalClass, 'intval')
    warning('validation is not rigorous without interval coefficient Scalars')
end

if nargin > 1
    printYZBounds = varargin{1}; % print Y, Z and error bounds to screen
else
    printYZBounds = false;
end

%% set up intval parameters, scaling, and truncation
mu1 = obj.Parameter(1); % these are the same numerical class as the chart
mu2 = obj.Parameter(2);
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
r1 = obj.Coordinate(5);
r2 = obj.Coordinate(5);

% generate base monomials of order 1
h11 = obj.Coordinate(1) + mu1; % x + mu1
h12 = obj.Coordinate(1) - mu2; % x - mu2
    
% generate base monomials of order 2 (no truncation)
h21 = fulltimes(h11, p); % (x + mu1)p
h22 = fulltimes(h12, p); % (x - mu2)p
yq = fulltimes(obj.Coordinate(3), obj.Coordinate(4)); % y*q
h23 = h21 + yq; % (x + mu1)p + yq
h24 = h22 + yq; % (x - mu2)p + yq
rr1 = fulltimes(obj.Coordinate(5), obj.Coordinate(5)); %r1^2
rr2 = fulltimes(obj.Coordinate(6), obj.Coordinate(6)); %r2^2

% generate base monomials of order 3 (no truncation)
rrr1 = fulltimes(obj.Coordinate(5), rr1); %r1^3
rrr2 = fulltimes(obj.Coordinate(5), rr2); %r1^3
h31 = 1-rrr1; % 1-r1^3
h32 = 1-rrr2; % 1-r2^3
    
% generate base monomials of order 4 (no truncation)
h41 = fulltimes(h11, h32); % (x + mu1)(1-r2^3)
h42 = fulltimes(h12, h31); % (x - mu2)(1-r1^3)
h43 = fulltimes(obj.Coordinate(3), h31); % y(1-r1^3)
h44 = fulltimes(obj.Coordinate(3), h32); % y(1-r2^3)

% generate base monomials of order 5 (no truncation)
h51 = fulltimes(rrr1, h24); % r1^3((x - mu2)p + yq)
h52 = fulltimes(rrr2, h23); % r2^3((x + mu1)p + yq)
    
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
Y2 = norm(tailProjection(int(2*q + h42*mu1 + h41*mu2)));
Y3 = norm(tailProjection(int(q)));
Y4 = norm(tailProjection(int(-2*p + h43*mu1 + h44*mu2)));
Y5 = norm(tailProjection(int(-1*h51)));
Y6 = norm(tailProjection(int(-1*h52)));
Y0 = Yscaling*max([Y1,Y2,Y3,Y4,Y5,Y6]) + obj.ErrorBound;

%% Z1 bounds (DERIVATIVE)
Zscaling = Yscaling/(M+1);
normVector = obj.Coordinate.norm(); % (||x||, ||p||, ||y||, ||q||, ||r1||, ||r2||) in R^6
nX = normVector(1);
nP = normVector(2);
nY = normVector(3);
nQ = normVector(4);
nR1 = normVector(5);
nR2 = normVector(6);

n1 = h11.norm(); % ||x + mu1||
n2 = h12.norm(); % ||x - mu2||
n3 = h31.norm(); % ||1 - r1^3||
n4 = h32.norm(); % ||1 - r2^3||

% Get derivative bounds
DT1 = intval([0,1,0,0,0,0]); 
DT2 = intval([mu1*n3 + mu2*n4, 0, 0, 2, 3*mu1*n2*nR1^2, 3*mu2*n1*nR2^2]); 
DT3 = intval([0,0,0,1,0,0]);
DT4 = intval([0, 2, mu1*n3 + mu2*n4, 3*mu1*nY*nR1^2, 3*mu2*nY*nR2^2]); 
DT5 = intval([nR1^3*nP, nR1^3*n2, nR1^3*nQ, nR1^3*nY, 3*nR1^2*h24.norm(), 0]);
DT6 = intval([nR2^3*nP, nR2^3*n1, nR2^3*nQ, nR2^3*nY, 0, 3*nR2^2*h23.norm()]);

ZMatrix = [DT1; DT2; DT3; DT4; DT5; DT6];
Z1 = Zscaling*norm(ZMatrix, Inf);

%% Z2 bounds (DERIVATIVE)
rStar = intval(10)*Y0;
nXStar = nX + rStar;
nPStar = nP + rStar;
nYStar = nY + rStar;
nQStar = nQ + rStar;
nR1Star = nR1 + rStar;
nR2Star = nR2 + rStar;
n1Star = nXStar + mu1; % mu1 + ||(x + rStar)||
n2Star = nXStar + mu2; % mu2 + ||(x + rStar)||
n3Star = 1 + nR1Star; % 1 + ||(r1 + rStar)^3||
n4Star = 1 + nR2Star; % 1 + ||(r2 + rStar)^3||
h23Star = n1Star*nPStar + nYStar*nQStar;
h24Star = n2Star*nPStar + nYStar*nQStar;


% Hessians of T_i evaluated at norm (||aBar|| + rStar)
H1 = intval(zeros(6,6)); 
H2 = intval([0, 0, 0, 0, 3*mu1*nR1Star^2, 3*mu2*nR2Star^2;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    3*mu1*nR1Star^2, 0, 0, 0, 6*mu1*nR1Star*n2Star, 0;...
    0, 0, 0, 0, 0, 0]); 
H3 = intval(zeros(6,6));
H4 = intval([0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 3*mu1*nR1Star^2, 3*mu2*nR2Star^2;...
    0, 0, 0, 0, 0, 0;...
    0, 0, 3*mu1*nR1Star^2, 0, 6*mu1*nR1Star*nYStar, 0;...
    0, 0, 3*mu2*nR2Star^2, 0, 0, 6*mu2*nR2Star*nYStar]); 
H5 = intval([0, nR1Star^3, 0, 0, 3*nR1Star^2*nPStar, 0;...
    nR1Star^3, 0, 0, 0, 3*nR1Star^2*n2Star, 0;...
    0, 0, 0, nR1Star^3, 3*nR1Star^2*nQStar, 0;...
    0, 0, nR1Star^3, 0, 3*nR1Star^2*nYStar, 0;...
    3*nR1Star^2*nPStar, 3*nR1Star^2*n2Star, 3*nR1Star^2*nQStar, 3*nR1Star^2*nYStar, 6*nR1Star*h24Star.norm(), 0;...
    0, 0, 0, 0, 0, 0]);
H6 = intval([0, nR2Star^3, 0, 0, 3*nR2Star^2*nPStar, 0;... 
    nR2Star^3, 0, 0, 0, 3*nR2Star^2*n1Star, 0;...
    0, 0, 0, nR2Star^3, 3*nR2Star^2*nQStar, 0;... 
    0, 0, nR2Star^3, 0, 3*nR2Star^2*nYStar, 0;...
    0, 0, 0, 0, 0, 0;...
    3*nR2Star^2*nPStar, 3*nR2Star^2*n1Star, 3*nR2Star^2*nQStar, 3*nR2Star^2*nYStar, 6*nR2Star*h23Star.norm(), 0]);

% Each element of ZVector is the norm of T(aBar + u) evaluated with the real scalar bound ||u|| = rStar. 
ZVector(1) = norm(H1, Inf); % ||D^2_T1(x,y)|| 
ZVector(2) = norm(H2, Inf); % ||D^2_T1(x,y)||
ZVector(3) = norm(H3, Inf); % ||D^2_T3(x,y)||
ZVector(4) = norm(H4, Inf); % ||D^2_T4(x,y)||
ZVector(5) = norm(H5, Inf); % ||D^2_T5(x,y)||
ZVector(6) = norm(H6, Inf); % ||D^2_T5(x,y)||
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
