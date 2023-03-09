function f = rk45regvectorfield(t,v,parameter,regType)
%RK45REGVECTORFIELD - Vector field evaluation for CRTBP field with regularization
%
%   Syntax:
%       f = RK45REGVECTORFIELD(t,v,mu,0) evaluation of f_0: R^6 ---> R^6 where mu is the small mass
%       f = RK45REGVECTORFIELD(t,v,[mu,C],0) evaluation of f_{1,2}: R^5 ---> R^5 where mu is the small mass and C is the regularization energy
%
%   Inputs:
%       parameter - MU or [MU,C] depending on regularization type
%       regType - 0, 1, or 2
%
%   Outputs:
%       f - evaluation of f_0,f_1, or f_2 at (t,x)
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 30-Mar-2019

with_automatic_diff = length(v) > 4;

switch regType
    case 0 % stanrdard CRTBP vector field
        if ~isequal(length(parameter), 1)
            error('The F0 field should only be passed a single parameter')
        end
        
        % unpack parameters
        MU = parameter;
        mu1 = 1 - MU;
        
        % unpack initial data
        x = v(1:6:end);
        p = v(2:6:end);
        y = v(3:6:end);
        q = v(4:6:end);
        
        if with_automatic_diff
            r1 = v(5:6:end);
            r2 = v(6:6:end);
        else  % added to allow evaluation of the vector field on R^4
            r1 = 1./sqrt((x - MU).^2 + y.^2);
            r2 = 1./sqrt((x + mu1).^2 + y.^2);
        end
        r1c = r1.^3;
        r2c = r2.^3;
        
        f = zeros(size(v));
        f1 = p;
        f2 = 2*q + mu1*(x-MU) + MU*(x+mu1) - mu1*(x-MU).*r1c - MU*(x+mu1).*r2c;
        f3 = q;
        f4 = -2*p + y - mu1*y.*r1c - MU*y.*r2c;
        f5 = -r1c.*((x-MU)*p + y.*q);
        f6 = -r2c.*((x+mu1)*p + y.*q);
        
        f(1:6:end) = f1;
        f(2:6:end) = f2;
        f(3:6:end) = f3;
        f(4:6:end) = f4;
        if with_automatic_diff
            f(5:6:end) = f5;
            f(6:6:end) = f6;
        end
        
    otherwise % regularized CRTBP vector field
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
        
        if ~isequal(length(parameter), 2)
            error('The F0 field should only be passed a single parameter')
        end
        
        % unpack parameters
        MU = parameter(1);
        C = parameter(2);
        
        if regType==1
            X_COORD = MU;
            OTHER_MASS = MU;
            REG_BIT = 1;
            
        elseif regType ==2
            X_COORD = MU-1;
            OTHER_MASS = 1-MU;
            REG_BIT = -1;
        end
        
        
        % unpack initial data
        x = v(1:5:end);
        p = v(2:5:end);
        y = v(3:5:end);
        q = v(4:5:end);
        if with_automatic_diff
            r = v(5:5:end);
        else  % added to allow evaluation of the vector field on R^4
            r = 1./sqrt((x.^2 + y.^2).^2 + 1 + 2*REG_BIT*(x.^2 - y.^2));
        end
        k = x.^2 + y.^2;
        
        f = zeros(size(v));
        f1 = p;
        f2 = 8*k.*q + 12*x.*k.^2 + 16*X_COORD*x.^3 + 4*(OTHER_MASS-C)*x + 8*OTHER_MASS*x.*(REG_BIT*x.^2 - 3*REG_BIT*y.^2 + 1).*r.^3;
        f3 = q;
        f4 =-8*k.*p + 12*y.*k.^2 - 16*X_COORD*y.^3 + 4*(OTHER_MASS-C)*y + 8*OTHER_MASS*y.*(-REG_BIT*y.^2 + 3*REG_BIT*x.^2 + 1).*r.^3;
        f5 = -2*r.^3.*(k.*(x.*p + y.*q) + REG_BIT*(x.*p - y.*q));
        
        f(1:5:end) = f1;
        f(2:5:end) = f2;
        f(3:5:end) = f3;
        f(4:5:end) = f4;
        if with_automatic_diff
            f(5:5:end) = f5;
        end
end

end % end rk45regvectorfield
