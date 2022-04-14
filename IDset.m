function initialData = IDset(mu, N, isValid)
%IDSET - Get a choice of initial data for testing CRTBP functions
%
%   Syntax:
%       initialData = IDSET(0,10,0) returns the line segment initial data set as degree 9 polynomials floating point coefficients.
%       initialData = IDSET(.25,20,0) returns the circle segment initial data set as degree 19 polynomials as floats.
%       initialData = IDSET(*,true) returns the testing initial data set as intervals
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 10-Apr-2019; Last revision: 19-May-2019

warning('This function is deprecated')
if isequal(mu,0)
    p1 = [0.7104, 0.2973, 0.2554, 0.2068];
    p2 = p1 + 1e-4*ones(size(p1));
    initialData = cat(2, 0.5*[p1+p2; p2-p1]', zeros(4,N-2));
    
else
    %  If mu <= 0.5 % set up some initial data on a circle in the zero energy level set of the F2 field otherwise the initial data lies on the F1 field
    thetaMid = pi/2;
    thetaRadius = pi;
    expCoeff = sqrt(8*mu)*exp(1i*thetaMid)*((1i*thetaRadius).^(0:N-1)./factorial(0:N-1));
    xInitial = zeros(1,N);
    pInitial = real(expCoeff);
    yInitial = zeros(1,N);
    qInitial = imag(expCoeff);
    initialData = [xInitial;pInitial;yInitial;qInitial];
end

if isValid
    initialData = intval(initialData);
end

end % end IDset

% Revision History:
%{
19 May 2019 - Added support for F1 field.
%}
