%BVPSOLVER - Refine a connection via a BVP formulation
%
%   Description:
%       BVPSOLVER description
%
%   Output:
%       BVPSOLVER output
%
%   Other m-files required: none
%   MAT-files required: none


%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 31-Dec-2019; Last revision: 31-Dec-2019

%% ================================================== SECTION 1 ==================================================

% INPUTS: 


%% ================================================== DEFINE A SINGLE NEWTON STEP OF THE BVP NEWTON METHOD ==================================================


% INPUTS: 
% alpha, beta, tau - real scalars
% t - vector of length n-1 which define fixed integration times for the multiple shooting segments
% U - 4-by-n array of column vectors, U = [u_1,...,u_{n-1}] which are initial points in phase space for the multiple shooting segments

n = size(U,2) + 1; % number of shooting segments



%% evalute g and Dg at (alpha, u_1,..., u_{n-1}, beta, tau)


% initialize evaluations
g = zeros(4, n); % g: R^{4n-1} --> R^{4n} is initialized as n-columns of component map evalutions in R^4
Dg = zeros(4*n, 4*n - 1); % Dg is a 4n-by-(4n-1) matrix

% construct flow variational time-tau maps
flowMap = @(tau, u)[SOME FUNCTION WITH MU]; % returns Phi(tau, u) in R^4 with mu fixed
variationalMap = @(tau, u)[SOME FUNCTION WITH MU]; % returns a vector in R^20. The first 4 components are Phi and the next 16 are 
% components of DPhi arranged in a column.


% construct first block of g and Dg
 
% iteratively construct middle blocks of g and Dg 
for j = 2:n-1 
    
end

% construct last block of g and Dg



 
%% ================================================== variationalsolver ==================================================

% INPUTS: 
% mu - small mass parameter
% x - length 20 vector. x = (x1, x2) where x1 in R^4 are the trajectory coordinates and x2 in R^16 is the variational 
%       solution shaped into a column vector


x1 = x(1:4); % usual cr3bp coordinates
Dphi = reshape(x(5:end),4,4); % reshape variational coordinates as a matrix
% Df = Dfcr3bp(mu, x1); 

% evaluate derivative of vector field
% define distance between (x,y) and the 2 primaries
P1 = [mu, 0]; % coordinates of large primary
P2 = [mu - 1, 0]; % coordinates of small primary
d1 = norm(x1([1,3]) - P1,2); % distance to large primary
d2 = norm(x1([1,3]) - P2,2); % distance to small primary

% define Hessian for Omega
d13 = d1.^(-3);
d15 = d1.^(-5);
d23 = d2.^(-3);
d25 = d2.^(-5);
H1 = 1 - (1 - mu)*(d13 - 3*(x1(1) - mu).^2.*d15) - mu*(d23 - 3*(x1(1) + 1 - mu).^2.*d25); % Omega_xx 
H2 = 3*x1(3).*((1 - mu)*(x1(1) - mu).*d15 + mu*(x1(1) + 1 - mu).*d25); % Omega_yx and Omega_xy
H3 = 1 - (1 - mu)*(d13 - 3*x1(3).^2.*d15) - mu*(d23 - 3*x1(3).^2.*d25); % Omega_yy
H = [H1, H2;H2, H3]; % Hessian of Omega

% define derivative
Df = zeros(4);
Df([5,15]) = 0;
Df(4,2) = -2;
Df(2,4) = 2;
Df([2,4],[2,4]) = H;
        
% evaluate vector field for variational equation and return as a column vector
F1 = cr3bp(mu, 1, x1); % evaluate usual vector field in first 4 coordinates
F2 = reshape(Df*Dphi, [], 1); % evaluate the variatonal equation as a column vector
F = [F1;F2]; % evaluation of extended vector field on R^20


 