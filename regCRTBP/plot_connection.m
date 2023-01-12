function plot_connection(connectionStruct, mu)
%PLOT_CONNECTION - Accepts a connection data structure (output by connection2BVP) and adds a plot of this connection in the configuration space of the current figure.
%
%   Syntax:
%       output = PLOT_CONNECTION(input)
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
%   Date: 12-Jan-2023;


for iStrand = 1:length(connectionStruct.Orbit)
    switch connectionStruct.RegVector(iStrand)
        case 0
            jSeg = connectionStruct.Orbit{iStrand};
            plot(jSeg(1,:), jSeg(3,:), 'g', 'LineWidth', 1)
            
        case 1
            jSeg = CRTBP2reg(connectionStruct.Orbit{iStrand}.', mu, -1).';
            plot(jSeg(1,:), jSeg(3,:), 'r', 'LineWidth', 1)
            
        case 2
            jSeg = CRTBP2reg(connectionStruct.Orbit{iStrand}.', mu, -2).';
            plot(jSeg(1,:), jSeg(3,:), 'b', 'LineWidth', 1)
    end
end
end % end plot_connection

