%CHECK_INVERSE - Verify the inverse function for Scalars
%
%   Description:
%       CHECK_INVERSE description
%
%   Output:
%       CHECK_INVERSE output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME
 
%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 31-Mar-2019; Last revision: 31-Mar-2019

a = Scalar(randi(10,5,5)-5, 'Taylor');
% a = Scalar(randi(10,5,1)-5,'Taylor');
a.Coefficient(1) = 3;
% da = a.dt;
% da.Coefficient
b = a.inv();
ab = a*b;
ab.Coefficient






% obj = a;


% ddtObj = obj.dt;
% c = ddtObj.Coefficient;
% a0 = obj.Coefficient(1);
% invObj = Scalar(1/a0, obj.Basis, obj.Truncation);
% for m = 1:obj.Truncation - 1
%     uu = invObj*invObj;
%     Aprime = c(m:-1:1);
%     newCoefficient = -(1/m)*dot(uu.Coefficient(1:m),Aprime);
%     invObj.Coefficient(m+1) = newCoefficient;
% end

% ab = obj*invObj;
% ab.Coefficient


% ab = a*b;
% ab.Coefficient

