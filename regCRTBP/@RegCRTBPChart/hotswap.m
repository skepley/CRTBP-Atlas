function newObj = hotswap(obj, swapTo)
%HOTSWAP - Changes between standard and regularized coordinates for a CRTBP boundary chart
%
%   HOTSWAP() - Maps a boundary chart in standard (F0) coordinates to a boundary chart in regularized (F1 or F2) coordinates or maps the inverse.
%               Hotswap is not defined to map directly between different regularized coordinates.
%
%   Syntax:
%       newBoundary = HOTSWAP(oldBoundary, k) % maps oldBd from Fi coordinates to a new boundary chart in Fk coordinates.
%       F2Boundary = HOTSWAP(HOTSWAP(F1Boundary, 0), 2) % Bootstraps the inverse mapping to map an F1 Boundary into an F2 Boundary.
%
%   Inputs:
%       obj - A boundary RegBRTBPChart in Fi coordinates
%       swapTo - k in {0,1,2}
%
%   Outputs:
%       newObj - A boundary RegBRTBPChart in Fk coordinates
%
%   Subfunctions: none
%   Classes required: RegCRTBPChart
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 14-Apr-2019; Last revision: 12-Jun-2019

disp('hotswapped')
obj.TimeSpan
mapDirection = swapTo - obj.RegType; % Possible values are {0, 1, -1, 2, -2}
mu = obj.Parameter(1);
coordinateBasis = obj.Coordinate(1).Basis;

switch mapDirection
    case 0 % Do not change coordinates i.e. this is the identity map.
        newObj = obj.deepcopy(); % Return a deep copy of the current boundary without changing any coordinates.
        return
    case 1 % map from F0 (standard) to F1 (regularized) coordinates
        % unpack the old coordinates
        x0 = obj.Coordinate(1);
        p0 = obj.Coordinate(2);
        y0 = obj.Coordinate(3);
        q0 = obj.Coordinate(4);
        C = CRTBPenergy([x0.Coefficient(1), p0.Coefficient(1), y0.Coefficient(1), q0.Coefficient(1)]', mu, 0);
        
        % map to regularized coordinates
        z = x0 + y0*1i;
        w = sqrt(z - mu);
        x1 = real(w);
        y1 = imag(w);
        p1 = 2*(x1*p0 + y1*q0);
        q1 = 2*(x1*q0 - y1*p0);
        newID = [x1.Coefficient; p1.Coefficient; y1.Coefficient; q1.Coefficient];
        newObj = RegCRTBPChart(newID, coordinateBasis, obj.TimeSpan, obj.Truncation, [mu,C], swapTo, 'InitialScaling', .1, 'boundary', true);
        
    case -1 % map from F1 (regularized) to F0 (standard) coordinates
        % unpack the old coordinates
        x1 = obj.Coordinate(1);
        p1 = obj.Coordinate(2);
        y1 = obj.Coordinate(3);
        q1 = obj.Coordinate(4);
        
        % set up the new coordinates
        xx = x1*x1;
        yy = y1*y1;
        k1 = xx + yy; % x^2 + y^2
        k2 = xx - yy; % x^2 - y^2
        k1Inverse = inv(k1); % compute 1/(x^2 + y^2)
        x0 = k2 + mu;
        p0 = (x1*p1 - y1*q1)*k1Inverse*0.5;
        y0 = 2*x1*y1;
        q0 = (y1*p1 + x1*q1)*k1Inverse*0.5;
        newID = [x0.Coefficient; p0.Coefficient; y0.Coefficient; q0.Coefficient];
        newObj = RegCRTBPChart(newID, coordinateBasis, obj.TimeSpan, obj.Truncation, mu, swapTo, 'InitialScaling', .1, 'boundary', true);
        
    case 2 % map from F0 (standard) to F2 (regularized) coordinates
        % unpack the old coordinates
        x0 = obj.Coordinate(1);
        p0 = obj.Coordinate(2);
        y0 = obj.Coordinate(3);
        q0 = obj.Coordinate(4);
        C = CRTBPenergy([x0.Coefficient(1), p0.Coefficient(1), y0.Coefficient(1), q0.Coefficient(1)]', mu, 0); % compute energy of regularization
        
        % map to regularized coordinates
        z = x0 + 1i*y0;
        w = sqrt(z - (mu-1));
        x2 = real(w);
        y2 = imag(w);
        p2 = 2*(x2*p0 + y2*q0);
        q2 = 2*(x2*q0 - y2*p0);
        newID = [x2.Coefficient; p2.Coefficient; y2.Coefficient; q2.Coefficient];
        newObj = RegCRTBPChart(newID, coordinateBasis, obj.TimeSpan, obj.Truncation, [mu,C], swapTo, 'InitialScaling', .1, 'boundary', true);
        
    case -2 % map from F2 (regularized) to F0 (standard) coordinates
        % unpack the old coordinates
        x1 = obj.Coordinate(1);
        p1 = obj.Coordinate(2);
        y1 = obj.Coordinate(3);
        q1 = obj.Coordinate(4);
        
        % set up the new coordinates
        xx = x1*x1;
        yy = y1*y1;
        k1 = xx + yy; % x^2 + y^2
        k2 = xx - yy; % x^2 - y^2
        invk1 = inv(k1); % compute 1/(x^2 + y^2)
        x0 = k2 + mu-1;
        p0 = (x1*p1 - y1*q1)*invk1*0.5;
        y0 = 2*x1*y1;
        q0 = (y1*p1 + x1*q1)*invk1*0.5;
        newID = [x0.Coefficient; p0.Coefficient; y0.Coefficient; q0.Coefficient];
        newObj = RegCRTBPChart(newID, coordinateBasis, obj.TimeSpan, obj.Truncation, mu, swapTo, 'InitialScaling', .1, 'boundary', true);
        
    otherwise
        error('Valid mapping values are {-2,-1,0,1,2}')
        
end % end mapDirection

newObj.Generation = obj.Generation; % update the new chart to the correct generation
newObj.ParentHandle = obj.ParentHandle; % update the new chart to the correct parent chart
newObj.SpatialSpan = obj.SpatialSpan; % update the new chart to the correct spatial domain
newObj.Crash = obj.Crash; % update the new chart to crash correctly


end % end hotswap

% Revision History:
%{
19 May 2019 - Added support for swapping from F0 to F1/F2.
12 Jun 2019 - Wrote documentation and added support for identity mapping to create deep copies.
%}
