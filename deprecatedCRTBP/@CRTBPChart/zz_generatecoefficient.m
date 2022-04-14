function generatecoefficient(obj, tau)
%GENERATECOEFFICIENT - recursively generate CRTBP system coefficients
%
%   Syntax:
%       obj.GENERATECOEFFICIENT() - compute Taylor coefficients and append as Scalar to obj.Coordinate
%
%   Description:
%       generatecoefficient(obj) - computer taylor series coefficients for the RegCRTBP vector field with regularized coordinates (x,p,y,q,r).
%           VECTOR FIELD:
%           xdot = p
%           pdot = 2q + mu1(x-mu2)(1-r1^3) + mu2(x+mu1)(1-r2^3)
%           ydot = q
%           qdot = -2p + mu1y(1-r1^3) + mu2y(1-r2^3)
%           r1dot = -r1^3((x-mu2)p + yq)
%           r2dot = -r2^3((x-mu2)p + yq)
%
%   Inputs:
%       obj - an instance of CRTBPChart
%       tau - specify vector field scaling
%
%   Classes required: @Chart, @Scalar, @CRTBPChart
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 20-Mar-2019; Last revision: 20-Mar-2019

if size(obj.Coordinate(1),1) > 1
    error('This chart has already been advected')
end

% retrieve parameters
mu1 = obj.Parameter(2);
mu2 = obj.Parameter(1);

% recursively compute finite approximation
for m = 1:obj.Truncation(1)-1
    scaleBy = tau/m;
    
    % generate base monomials of order 1
    h11 = obj.Coordinate(1) + mu1; % x + mu1
    h12 = obj.Coordinate(1) - mu2; % x - mu2
    
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
    rrr2 = obj.Coordinate(5)*rr2; %r1^3
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
            % evaluate vector field as matlab array
            %             xDotNext = scaleBy*
            %             pDotNext = scaleBy*
            %             yDotNext = scaleBy*
            %             qDotNext = scaleBy*
            %             r1DotNext = scaleBy*
            %             r2DotNext = scaleBy*
        case 2
            % evaluate vector field as matlab array
            xDotNext = scaleBy*(obj.Coordinate(2).Coefficient(m,:));
            pDotNext = scaleBy*(2*obj.Coordinate(4).Coefficient(m,:) + mu1*h42.Coefficient(m,:) + mu2*h41.Coefficient(m,:));
            yDotNext = scaleBy*(obj.Coordinate(4).Coefficient(m,:));
            qDotNext = scaleBy*(-2*obj.Coordinate(2).Coefficient(m,:) + mu1*h43.Coefficient(m,:) + mu2*h44.Coefficient(m,:));
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
    
end % for
end % generatecoefficient


% Revision History:
%{

%}
