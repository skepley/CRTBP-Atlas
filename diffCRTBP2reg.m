function DFInverse = diffCRTBP2reg(data)
%CRTBP2REG - Evaluate differential of the inverse regularization maps which take regularized (F1/F2) CRTBP coordinates
% to standard (F0) CRTBP coordinates.
%
%   Syntax:
%       newCoordinate = diffCRTBP2REG(data, mu) returns the derivative of the F_i ---> F_0 coordinate map.
%
%   Inputs:
%       data - A single vector in R^4 for coordinates u_i = (xi,pi,yi,qi) for i = 1, or 2. This function is NOT vectorized. 
%
%   Outputs:
%       DFInverse - The derivative (as a 4-by-by matrix) of the inverse map taking u_i |--> u_0.
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 10-Mar-2021; Last revision: 10-Mar-2021

mapDirection = -1;  % remove this if the forward derivative is implemented
if isequal(mapDirection, -1)  % compute derivative of inverse regulatization map
    % unpack coordinates (in X1 or X2) and set parameters
    x = data(1);
    p = data(2);
    y = data(3);
    q = data(4);
    
    % Construct Jacobian for F^-1(x,p,y,q)
    D = x.^2 + y.^2;
    DFInverse = zeros(4);
    DFInverse(1, :) = [2*x, 0, -2*y, 0];
    DFInverse(2,:) = [(p*(y.^2 - x.^2) + 2*x*y*q)/(2*D.^2), x./(2*D), (q*(y.^2 - x.^2) - 2*x*y*p)/(2*D.^2), -y./(2*D)];
    DFInverse(3,:) = [2*y, 0, 2*x, 0]; 
    DFInverse(4,:) = [(q*(x.^2 -y.^2) - 2*x*y*p)/(2*D.^2), y./(2*D), (p*(x.^2 - y.^2) + 2*x*y*q)/(2*D.^2), -x./(2*D)];


else
    error('Derivative of forward regularization map not yet implemented')
end


end % end diffCRTBP2reg

% Revision History:
%{

%}
