function [output1,output2] = equalmassidealdomain(input1,input2,input3,varargin)
%EQUALMASSIDEALDOMAIN - One line description of what the function or script performs (H1 line)
%
%   EQUALMASSIDEALDOMAIN() - A more detailed description of the function
%
%   Syntax:
%       output = EQUALMASSIDEALDOMAIN(input1, input2)
%       [output1, output2] = EQUALMASSIDEALDOMAIN(input1, input2, input3)
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
%   Date: 17-Apr-2020; Last revision: 17-Apr-2020

%% parse input 
% p = inputParser;
% addRequired(p, 'input1')
% addRequired(p, 'input2')
% addRequired(p, 'input3')
% addParameter(p, 'Parameter1', default1)
% addParameter(p, 'Parameter2', default2)
% addParameter(p, 'Parameter3', default3)
% 
% parse(p,input1,input2,input3,varargin{:})
% parameter1 = p.Results.Parameter1;
% parameter2 = p.Results.Parameter2;
% parameter3 = p.Results.Parameter3;


load ideal_domain_boundary_C3_equalmass
mu = 1/2;


end % end equalmassidealdomain

% Revision History:
%{

%}
