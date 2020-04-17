function generatecoefficient(obj, tau)
%GENERATECOEFFICIENT - recursively generate RegCRTBP system coefficients
%
%   Syntax:
%       obj.GENERATECOEFFICIENT() - compute Taylor coefficients and append as Scalar to obj.Coordinate
%
%   Description:
%       generatecoefficient(obj) - computer taylor series coefficients for the RegCRTBP vector field with standard or regularized coordinates.
%           F_1 OR F_2 VECTOR FIELD WITH VARIABLES (x,p,y,q,r)
%           xdot = p
%           pdot = 8(x^2 + y^2)q + 12x(x^2 + y^2)^2 + 16[X_COORD]x^3 + 4([OTHER_MASS] - C)x + 8[OTHER_MASS]x([REG_BIT]x^2 - 3[REG_BIT]y^2 + 1)r^3
%           ydot = q
%           qdot = 8(x^2 + y^2)p + 12y(x^2 + y^2)^2 - 16[X_COORD]y^3 + 4([OTHER_MASS] - C)y + 8[OTHER_MASS]y(-[REG_BIT]y^2 + 3[REG_BIT]x^2 + 1)r^3
%           rdot = -2r^3((x^2 + y^2)(xp + yq) +[REG_BIT](xp - yq))
%           where
%           [X_COORD] is the x-coordinate of the primary being desingularized (so either mu or mu-1)
%           [OTHER_MASS] is the mass of the opposite primary which is not being desingularized (so either 1-mu or mu)
%           [REG_BIT] = 1 (desingularize the large primary) or -1 (desingularize the small primary)
%
%           F_0 VECTOR FIELD: WITH VARIABLES (x,p,y,q,R1,R2)
%           xdot = p
%           pdot = 2q + mu1(x-mu2)(1-r1^3) + mu2(x+mu1)(1-r2^3)
%           ydot = q
%           qdot = -2p + mu1y(1-r1^3) + mu2y(1-r2^3)
%           r1dot = -r1^3((x-mu2)p + yq)
%           r2dot = -r2^3((x+mu1)p + yq)
%           where
%           mu2 = mu is the small mass and mu1 = 1-mu is the large mass
%
%   Inputs:
%       obj - an instance of RegCRTBPChart
%       tau - specify vector field scaling
%
%   Classes required: @Chart, @Scalar, @RegCRTBPChart
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 07-Mar-2019; Last revision: 20-Feb-2020

if size(obj.Coordinate(1),1) > 1
    error('This chart has already been advected')
end

switch obj.RegType
    case 0 % recursively compute Taylor coefficients for standard CRTBP vector field
        
        % retrieve parameters
        MU = obj.Parameter(1); % small mass
        MU1 = 1-MU; % large mass
        
        % recursively compute finite approximation
        for m = 1:obj.Truncation(1)-1
            scaleBy = tau/m;
            
            % generate base monomials of order 1
            h11 = obj.Coordinate(1) + MU1; % x + mu1
            h12 = obj.Coordinate(1) - MU; % x - mu2
            
            % generate base monomials of order 2
            h21 = h11*obj.Coordinate(2); % (x + mu1)p
            h22 = h12*obj.Coordinate(2); % (x - mu2)p
            yq = obj.Coordinate(3)*obj.Coordinate(4); % y*q
            h23 = h21 + yq; % (x + mu1)p + yq
            h24 = h22 + yq; % (x - mu2)p + yq
            rr1 = obj.Coordinate(5)*obj.Coordinate(5); %r1^2
            rr2 = obj.Coordinate(6)*obj.Coordinate(6); %r2^2
            
            % generate base monomials of order 3
            rrr1 = obj.Coordinate(5)*rr1; %r1^3
            rrr2 = obj.Coordinate(6)*rr2; %r2^3
            h31 = 1-rrr1; % 1-r1^3
            h32 = 1-rrr2; % 1-r2^3
            
            % generate base monomials of order 4
            h41 = h11*h32; % (x + mu1)(1-r2^3)
            h42 = h12*h31; % (x - mu2)(1-r1^3)
            h43 = obj.Coordinate(3)*h31; % y(1-r1^3)
            h44 = obj.Coordinate(3)*h32; % y(1-r2^3)
            
            % generate base monomials of order 5
            h51 = rrr1*h24; % r1^3((x - mu2)p + yq)
            h52 = rrr2*h23; % r2^3((x + mu1)p + yq)
            
            switch obj.Dimension(1) % Surface dimension
                case 1
                    xDotNext = scaleBy*(obj.Coordinate(2).Coefficient(m));
                    pDotNext = scaleBy*(2*obj.Coordinate(4).Coefficient(m) + h42.Coefficient(m)*MU1 + h41.Coefficient(m)*MU);
                    yDotNext = scaleBy*(obj.Coordinate(4).Coefficient(m));
                    qDotNext = scaleBy*(-2*obj.Coordinate(2).Coefficient(m) + h43.Coefficient(m)*MU1 + h44.Coefficient(m)*MU);
                    r1DotNext = scaleBy*(-h51.Coefficient(m));
                    r2DotNext = scaleBy*(-h52.Coefficient(m));
                    
                case 2 % initial data is 1-dimensional
                    xDotNext = scaleBy*(obj.Coordinate(2).Coefficient(m,:));
                    pDotNext = scaleBy*(2*obj.Coordinate(4).Coefficient(m,:) + h42.Coefficient(m,:)*MU1 + h41.Coefficient(m,:)*MU);
                    yDotNext = scaleBy*(obj.Coordinate(4).Coefficient(m,:));
                    qDotNext = scaleBy*(-2*obj.Coordinate(2).Coefficient(m,:) + h43.Coefficient(m,:)*MU1 + h44.Coefficient(m,:)*MU);
                    r1DotNext = scaleBy*(-h51.Coefficient(m,:));
                    r2DotNext = scaleBy*(-h52.Coefficient(m,:));
                otherwise
                    error('Not implemented')
            end
            % append next coefficient
            obj.Coordinate(1) = append(obj.Coordinate(1), xDotNext);
            obj.Coordinate(2) = append(obj.Coordinate(2), pDotNext);
            obj.Coordinate(3) = append(obj.Coordinate(3), yDotNext);
            obj.Coordinate(4) = append(obj.Coordinate(4), qDotNext);
            obj.Coordinate(5) = append(obj.Coordinate(5), r1DotNext);
            obj.Coordinate(6) = append(obj.Coordinate(6), r2DotNext);
        end % for loop
        
    otherwise % recursively compute Taylor coefficients for regularized CRTBP vector field
        
        % retrieve parameters
        MU = obj.Parameter(1);
        C = obj.Parameter(2);
        if isequal(obj.RegType,1) % desingularize large primary
            X_COORD = MU;
            OTHER_MASS = MU;
            REG_BIT = 1;
            
        else % desingularize small primary
            X_COORD = MU-1;
            OTHER_MASS = 1-MU;
            REG_BIT = -1;
        end
        
        % recursively compute finite approximation
        for m = 1:obj.Truncation(1)-1
            scaleBy = tau/m;
            
            % generate base monomials of order 2
            xx = obj.Coordinate(1)*obj.Coordinate(1); % x^2
            yy = obj.Coordinate(3)*obj.Coordinate(3); % y^2
            k1 = xx + yy; % x^2 + y^2
            rr = obj.Coordinate(5)*obj.Coordinate(5); % r^2
            xp = obj.Coordinate(1)*obj.Coordinate(2); % x*p
            yq = obj.Coordinate(3)*obj.Coordinate(4); % y*q
            h1 = xp + yq; % x*p + y*q
            h2 = xp - yq; % x*p - y*q
            g21 = REG_BIT*xx - 3*REG_BIT*yy + 1; % [REG_BIT]x^2 - 3[REG_BIT]y^2 + 1
            g22 = 3*REG_BIT*xx - REG_BIT*yy + 1; % -[REG_BIT]y^2 + 3[REG_BIT]x^2 + 1
            
            % generate base monomials of order 3
            xxx = obj.Coordinate(1)*xx; % x^3
            yyy = obj.Coordinate(3)*yy; % y^3
            rrr = obj.Coordinate(5)*rr; % r^3
            k1q = k1*obj.Coordinate(4); % (x^2 + y^2)q
            k1p = k1*obj.Coordinate(2); % (x^2 + y^2)p
            
            % generate base monomials of order 4
            xrrr = obj.Coordinate(1)*rrr; % xr^3
            yrrr = obj.Coordinate(3)*rrr; % yr^3
            k1k1 = k1*k1; % (x^2 + y^2)^2
            k1h1 = k1*h1; % (x^2 + y^2)(xp + pq)
            g41 = k1h1 + REG_BIT*h2; % (x^2 + y^2)(xp + pq) + [REG_BIT]*(xp - pq)
            
            % generate base monomials of order 5
            xk1k1 = k1k1*obj.Coordinate(1); % x(x^2 + y^2)^2
            yk1k1 = k1k1*obj.Coordinate(3); % y(x^2 + y^2)^2
            
            % generate base monomials of order 6
            g61 = g21*xrrr; % (x^2 - 3y^2 + 1)xr^3
            g62 = g22*yrrr; % (-y^2 + 3x^2 + 1)yr^3
            
            % generate base monomials of order 7
            g71 = rrr*g41; %  r^3((x^2 + y^2)(xp + pq) - (xp - pq))
            
            switch obj.Dimension(1) % Surface dimension
                case 1
                    % evaluate vector field as matlab array
                    xDotNext = scaleBy*(obj.Coordinate(2).Coefficient(m)); % xDot = p
                    pDotNext = scaleBy*(8*k1q.Coefficient(m) + 12*xk1k1.Coefficient(m) + 16*X_COORD*xxx.Coefficient(m) +...
                        + 4*(OTHER_MASS-C)*obj.Coordinate(1).Coefficient(m) + 8*OTHER_MASS*g61.Coefficient(m));
                    yDotNext = scaleBy*(obj.Coordinate(4).Coefficient(m)); % yDot = q
                    qDotNext = scaleBy*(-8*k1p.Coefficient(m) + 12*yk1k1.Coefficient(m) - 16*X_COORD*yyy.Coefficient(m) +...
                        + 4*(OTHER_MASS-C)*obj.Coordinate(3).Coefficient(m) + 8*OTHER_MASS*g62.Coefficient(m));
                    rDotNext = scaleBy*(-2*g71.Coefficient(m));
                    
                case 2
                    % evaluate vector field as matlab array
                    xDotNext = scaleBy*(obj.Coordinate(2).Coefficient(m,:)); % xDot = p
                    pDotNext = scaleBy*(8*k1q.Coefficient(m,:) + 12*xk1k1.Coefficient(m,:) + 16*X_COORD*xxx.Coefficient(m,:) +...
                        + 4*(OTHER_MASS-C)*obj.Coordinate(1).Coefficient(m,:) + 8*OTHER_MASS*g61.Coefficient(m,:));
                    yDotNext = scaleBy*(obj.Coordinate(4).Coefficient(m,:)); % yDot = q
                    qDotNext = scaleBy*(-8*k1p.Coefficient(m,:) + 12*yk1k1.Coefficient(m,:) - 16*X_COORD*yyy.Coefficient(m,:) +...
                        + 4*(OTHER_MASS-C)*obj.Coordinate(3).Coefficient(m,:) + 8*OTHER_MASS*g62.Coefficient(m,:));
                    rDotNext = scaleBy*(-2*g71.Coefficient(m,:));
                otherwise
                    error('Not implemented')
            end
            
            % append next coefficient
            obj.Coordinate(1) = append(obj.Coordinate(1), xDotNext);
            obj.Coordinate(2) = append(obj.Coordinate(2), pDotNext);
            obj.Coordinate(3) = append(obj.Coordinate(3), yDotNext);
            obj.Coordinate(4) = append(obj.Coordinate(4), qDotNext);
            obj.Coordinate(5) = append(obj.Coordinate(5), rDotNext);
            
        end % for
        
end % switch
end % generatecoefficient


% Revision History:
%{
09-Apr-2019 - Merged all 3 types of CRTBP chart into a single class. This file computes coefficients for any of these 3.
20-Feb-2010 - Added support for 0 dimensional initial data.
%}
