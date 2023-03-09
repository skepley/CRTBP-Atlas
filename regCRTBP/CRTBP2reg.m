function newCoordinate = CRTBP2reg(data, mu, mapDirection)
%CRTBP2REG - Change of coordinates between physical (F0) CRTBP coordinates and regularized (F1/F2) CRTBP coordinates
%
%   Syntax:
%       newCoordinate = CRTBP2REG(data, mu, 1) returns the F1 regularized coordinates for data specified in F0 coordinates
%       newCoordinate = CRTBP2REG(data, mu, -1) returns the F0 coordinates for data specified in F1 coordinates
%       newCoordinate = CRTBP2REG(data, mu, 2) returns the F2 regularized coordinates for data specified in F0 coordinates
%       newCoordinate = CRTBP2REG(data, mu, -2) returns the F0 coordinates for data specified in F2 coordinates
%
%   Inputs:
%       mapDirection: -1, 1, -2, or 2 to compute h1 or h2 conjugacy maps or their inverses
%       mu - The mass of the small primary satisfying 0 < mu <= 0.5
%       data - An M-by-4 array of row vectors in R^4 for coordinates (xi,pi,yi,qi) for i = 0,1, or 2. 
%
%   Outputs:
%       newCoordinate - The image of the point in the new coordinates in R^6 or an array of points as an M-by-6 array
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 04-Apr-2019; Last revision: 12-June-2019

mu1 = 1-mu; % large mass value

% force shape of data to be specified as row vectors of points in R^4
if ~isequal(size(data,2),4)
    error('Data must have 4 columns (as a matrix of row vectors specifying points in R^4)')
end

switch mapDirection
    case 1 % map from F0 (physical) to F1 (regularized) coordinates
        % unpack coordinates
        x0 = data(:,1); 
        p0 = data(:,2);
        y0 = data(:,3);
        q0 = data(:,4);
        
        % map to regularized coordinates
        z = x0 + 1i*y0;
        w = sqrt(z - mu);
        x1 = real(w);
        y1 = imag(w);
        p1 = 2*(x1.*p0 + y1.*q0);
        q1 = 2*(x1.*q0 - y1.*p0);
        
        % map automatic differentation coordinates and energy
        k1 = x1.^2 + y1.^2;
        k2 = x1.^2 - y1.^2;
        r1 = 1./sqrt(k1.^2 + 1 + 2*k2);
        C = CRTBPenergy(data(1,:), mu, 0); % compute energy of regularization
        newCoordinate = [x1,p1,y1,q1,r1,C*ones(size(x1))];
        
    case -1 % map from F1 (regularized) to F0 (physical) coordinates
        % unpack coordinates
        x1 = data(:,1);
        p1 = data(:,2);
        y1 = data(:,3);
        q1 = data(:,4);
        
        % map to CRTBP coordinates
        N = x1.^2 + y1.^2;
        x0 = x1.^2 - y1.^2 + mu;
        p0 = (x1.*p1 - y1.*q1)./(2*N);
        y0 = 2*x1.*y1;
        q0 = (y1.*p1 + x1.*q1)./(2*N);
        
        % map automatic differentation coordinates
        r0 = 1./sqrt((x0 - mu).^2 + y0.^2);
        s0 = 1./sqrt((x0 + mu1).^2 + y0.^2); 
        newCoordinate = [x0,p0,y0,q0,r0,s0];
        
    case 2 % map from F0 (physical) to F2 (regularized) coordinates
        % unpack coordinates
        x0 = data(:,1); 
        p0 = data(:,2);
        y0 = data(:,3);
        q0 = data(:,4);
        
        % map to regularized coordinates
        z = x0 + 1i*y0;
        w = sqrt(z - (mu-1));
        x2 = real(w);
        y2 = imag(w);
        p2 = 2*(x2.*p0 + y2.*q0);
        q2 = 2*(x2.*q0 - y2.*p0);
        
        % map automatic differentation coordinates and energy
        k1 = x2.^2 + y2.^2;
        k2 = x2.^2 - y2.^2;
        r2 = 1./sqrt(k1.^2 + 1 + 2*k2);
        C = CRTBPenergy(data(1,:), mu, 0); % compute energy of regularization
        newCoordinate = [x2,p2,y2,q2,r2,C*ones(size(x2))];
        
    case -2 % map from F2 (regularized) to F0 (physical) coordinates
        % unpack coordinates
        x2 = data(:,1);
        p2 = data(:,2);
        y2 = data(:,3);
        q2 = data(:,4);
        
        % map to CRTBP coordinates
        N = x2.^2 + y2.^2;
        x0 = x2.^2 - y2.^2 + mu - 1;
        p0 = (x2.*p2 - y2.*q2)./(2*N);
        y0 = 2*x2.*y2;
        q0 = (y2.*p2 + x2.*q2)./(2*N);
        
        % map automatic differentation coordinates
        r0 = 1./sqrt((x0 - mu).^2 + y0.^2);
        s0 = 1./sqrt((x0 + mu1).^2 + y0.^2);
        newCoordinate = [x0,p0,y0,q0,r0,s0];
        
    otherwise
        error('invalid map direction')
end
end % end CRTBP2reg

% Revision History:
%{
19 May 2019 - Fixed this function to perform any coordinate changes and updated documentation accordingly. 
12 Jun 2019 - Fixed vectorization and added support for intvals to be used with ruling out intersections via ell-one boxes.
%}
