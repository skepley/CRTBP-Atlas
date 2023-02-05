function [output1,output2] = verify_true_orbit(rawConnectionData, parsedConnectionData, mu)
%VERIFY_TRUE_ORBIT - plot a connection over top of its continuation through the intersection to verify visually whether it is 
% a true orbit or pseudo orbit. 
%
%   Syntax:
%       output = VERIFY_TRUE_ORBIT(input)
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
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
%   email: s.kepley@vu.nl
%   Date: 24-Jan-2023; 

[sData, uData] = continue_through_intersection(rawConnectionData);
figure 
hold on
% plot_primaries(mu, L4)
plot(sData(:, 1), sData(:, 2), 'r')
plot(uData(:, 1), uData(:, 2), 'k')
plot_connection(parsedConnectionData, mu);
end % end verify_true_orbit
