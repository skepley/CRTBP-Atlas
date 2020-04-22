%CHECK_TAYLORREGTIME - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       CHECK_TAYLORREGTIME description
%
%   Output:
%       CHECK_TAYLORREGTIME output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 22-Apr-2020; Last revision: 22-Apr-2020

%% ================================================== SECTION 1 ==================================================
p = [3,6,-4,2]; % 3*t^3 + 6*t^2 - 4*t + 2
pp = conv(p,p);
ppI = polyint(pp); % anitderivative of p^2
ground_truth = diff(polyval(ppI, [1,0]))

% my method
d = size(p, 2);
chk = sum(conv(p,p)./(2*d-1:-1:1))




%% ================================================== SECTION 2 ==================================================




