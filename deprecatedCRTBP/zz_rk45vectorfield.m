function f = rk45vectorfield(t,v,mu)
%RK45VECTORFIELD - One line description of what the function or script performs (H1 line)
%
%   RK45VECTORFIELD() - A more detailed description of the function
%
%   Syntax:
%       output = RK45VECTORFIELD(input1, input2)
%       [output1, output2] = RK45VECTORFIELD(input1, input2, input3)
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Apr-2019; Last revision: 04-Apr-2019

% unpack initial data
x = v(1:6:end);
p = v(2:6:end);
y = v(3:6:end);
q = v(4:6:end);
r1 = v(5:6:end);
r2 = v(6:6:end);
mu1 = 1 - mu;

r1c = r1.^3;
r2c = r2.^3;

f = zeros(size(v));
f1 = p;
f2 = 2*q + mu1*(x-mu) + mu*(x+mu1) - mu1*(x-mu).*r1c - mu*(x+mu1).*r2c;
f3 = q;
f4 = -2*p + y - mu1*y.*r1c - mu*y.*r2c;
f5 = -r1c.*((x-mu)*p + y.*q);
f6 = -r2c.*((x+mu1)*p + y.*q);

f(1:6:end) = f1;
f(2:6:end) = f2;
f(3:6:end) = f3;
f(4:6:end) = f4;
f(5:6:end) = f5;
f(6:6:end) = f6;

end % end rk45vectorfield

% Revision History:
%{

%}
