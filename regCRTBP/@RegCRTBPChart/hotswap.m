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
%   Date: 14-Apr-2019; Last revision: 17-Jun-2020

obj.TimeSpan
mu = obj.Parameter(1);
coordinateBasis = obj.Coordinate(1).Basis;

if isequal(swapTo, 0) || isequal(obj.RegType, 0)  % this is the original function which mapped only between standard and regularized coordinates
    mapDirection = swapTo - obj.RegType; % Possible values are {0, 1, -1, 2, -2}
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
            %             C = CRTBPenergy([x0.Coefficient(1), p0.Coefficient(1), y0.Coefficient(1), q0.Coefficient(1)]', mu, 0);
            C = meanenergy(obj);  % compute mean anergy evaluated over 100 uniformly spaced nodes.

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
            %             C = CRTBPenergy([x0.Coefficient(1), p0.Coefficient(1), y0.Coefficient(1), q0.Coefficient(1)]', mu, 0); % compute energy of regularization at center of boundary arc
            C = meanenergy(obj);  % compute mean anergy evaluated over 100 uniformly spaced nodes.
            
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
else  % updated version of the code includes swapping between regularized coordinates directly
    switch obj.RegType
        case 1  % compute u_2 = F_2 o F_1^{-1} (u_1)
            % unpack the old coordinates
            x1 = obj.Coordinate(1);
            p1 = obj.Coordinate(2);
            y1 = obj.Coordinate(3);
            q1 = obj.Coordinate(4);
            C = obj.Parameter(2);
            
            % compute w and z
            xx = x1*x1;
            yy = y1*y1;
            z2 = xx - yy + 1 + 1i*2*x1*y1;  % x0+mu + i*y0 expressed in terms of x1,y1
            w2 = sqrt(z2);
            x2 = real(w2);
            y2 = imag(w2);
            
            % compute products and sums only one time
            xp = x1*p1;
            yq = y1*q1;
            xq = x1*q1;
            yp = y1*p1;
            xp_yq = xp - yq;
            yp_xq = yp + xq;
            k1 = xx + yy;
            invk1 = inv(k1);
            p0Term = xp_yq*invk1;  %(x1p1 - y1q1)/(x1^2 + y1^2)
            q0Term = yp_xq*invk1;  %(y1p1 + x1q1)/(x1^2 + y1^2)
            
            % use previous computations x2,y2, and products/sums to compute p2,q2
            p2 = p0Term*x2 + q0Term*y2;
            q2 = q0Term*x2 - p0Term*y2;
            newID = [x2.Coefficient; p2.Coefficient; y2.Coefficient; q2.Coefficient];
            newObj = RegCRTBPChart(newID, coordinateBasis, obj.TimeSpan, obj.Truncation, [mu,C], swapTo, 'InitialScaling', .1, 'boundary', true);
            
        case 2  % compute u_1 = F_1 o F_2^{-1} (u_2)
            % unpack the old coordinates
            x2 = obj.Coordinate(1);
            p2 = obj.Coordinate(2);
            y2 = obj.Coordinate(3);
            q2 = obj.Coordinate(4);
            C = obj.Parameter(2);
            
            
            % compute w and z
            xx = x2*x2;
            yy = y2*y2;
            z1 = xx - yy - 1 + 1i*2*x2*y2;  % x0+mu-1 + i*y0 expressed in terms of x2,y2
            w1 = sqrt(z1);
            x1 = real(w1);
            y1 = imag(w1);
            
            % compute products and sums only one time
            xp = x2*p2;
            yq = y2*q2;
            xq = x2*q2;
            yp = y2*p2;
            xp_yq = xp - yq;
            yp_xq = yp + xq;
            k1 = xx + yy;
            invk1 = inv(k1);
            p0Term = xp_yq*invk1;  %(x2p2 - y2q2)/(x2^2 + y2^2)
            q0Term = yp_xq*invk1;  %(y2p2 + x2q2)/(x2^2 + y2^2)
            
            % use previous computations x1,y1, and products/sums to compute p1,q1
            p1 = p0Term*x1 + q0Term*y1;
            q1 = q0Term*x1 - p0Term*y1;
            newID = [x1.Coefficient; p1.Coefficient; y1.Coefficient; q1.Coefficient];
            newObj = RegCRTBPChart(newID, coordinateBasis, obj.TimeSpan, obj.Truncation, [mu,C], swapTo, 'InitialScaling', .1, 'boundary', true);
            
    end
    
end

newObj.Generation = obj.Generation; % update the new chart to the correct generation
newObj.ParentHandle = obj.ParentHandle; % update the new chart to the correct parent chart
newObj.SpatialSpan = obj.SpatialSpan; % update the new chart to the correct spatial domain
newObj.Crash = obj.Crash; % update the new chart to crash correctly


end % end hotswap


function meanEnergy = meanenergy(chart)
%MEANENERGY - compute the Mean energy along an F0 boundary chart discretzied at 100 uniform points

s = linspace(-1,1,100).';
chartEval = cell2mat(chart.eval(s));
energy = CRTBPenergy(chartEval(:, 1:4), chart.Parameter, chart.RegType);
meanEnergy = mean(energy);
end


% Revision History:
%{
19 May 2019 - Added support for swapping from F0 to F1/F2.
12 Jun 2019 - Wrote documentation and added support for identity mapping to create deep copies.
17 Jun 2020 - Added support for coordinate transform from F1 to F2 and F2 to F1. Also changed the computation of the regularization energy
when mapping from F0 to F1/F2 to take the mean energy across the boundary arc.
%}
