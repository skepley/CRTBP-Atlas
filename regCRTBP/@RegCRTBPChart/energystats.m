function [meanData, stdData] = energystats(obj, nPts, varargin)
%ENERGYSTATS - Compute mean and std of the energy of points lying in this chart
%
%   ENERGYSTATS() - A more detailed description of the function
%
%   Syntax:
%       output = ENERGYSTATS(input1, input2)
%       [output1, output2] = ENERGYSTATS(input1, input2, input3)
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
%   email: s.kepley@vu.nl
%   Date: 08-Jul-2019; Last revision: 08-Jul-2019

%% parse input 
s = linspace(-1,1,nPts);
t = linspace(0,1,nPts);
[S,T] = meshgrid(s,t);
data = [reshape(S, [], 1), reshape(T, [], 1)];
chartEval = cell2mat(obj.eval(data));
energy = CRTBPenergy(chartEval(:, 1:4), obj.Parameter, obj.RegType);
meanData = mean(energy);
stdData = std(energy);
% statData = [mean(energy), std(energy)];

end % end energystats

% Revision History:
%{

%}
