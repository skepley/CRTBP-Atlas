function [output1,output2] = plot(obj, nNode, varargin)
%PLOT - One line description of what the function or script performs (H1 line)
%
%   Syntax:
%       output = PLOT(input)
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
%   Date: 23-Feb-2023; 

sampleData = obj.sample(nNode);
mu = obj.Parameter;
hold_state = ishold;  % check if hold is on 
hold on
for iStrand = 1:length(sampleData.Orbit)
    switch sampleData.RegVector(iStrand)
        case 0
            jSeg = sampleData.Orbit{iStrand};
            plot(jSeg(1,:), jSeg(3,:), 'g', 'LineWidth', 1, varargin{:})
            
        case 1
            jSeg = CRTBP2reg(sampleData.Orbit{iStrand}.', mu, -1).';
            plot(jSeg(1,:), jSeg(3,:), 'r', 'LineWidth', 1, varargin{:})
            
        case 2
            jSeg = CRTBP2reg(sampleData.Orbit{iStrand}.', mu, -2).';
            plot(jSeg(1,:), jSeg(3,:), 'b', 'LineWidth', 1, varargin{:})
    end
end
if ~hold_state
    hold off  % turn hold back off if it was off when we started
end

end % end plot

