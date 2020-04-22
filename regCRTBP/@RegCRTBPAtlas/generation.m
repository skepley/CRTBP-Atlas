function thisGeneration = generation(obj, generationIdx)
%GENERATION - Return the cohort charts which lie in the generation index
%
%   GENERATION() - A more detailed description of the function
%
%   Syntax:
%       output = GENERATION(input1, input2)
%       [output1, output2] = GENERATION(input1, input2, input3)
%
%   Inputs:
%       generationIdx - A vector of integers describing which generations to return
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
%   Date: 06-Apr-2019; Last revision: 06-Apr-2019


fullGenerationIdx = [obj.Chart.Generation]; % list of generation indices for entire atlas
thisGeneration = obj.Chart(ismember(fullGenerationIdx,generationIdx)); % return the charts in the generationIdx specified.

end % end generation

% Revision History:
%{

%}
