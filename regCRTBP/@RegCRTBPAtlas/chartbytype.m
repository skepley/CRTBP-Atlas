function chartSlice = chartbytype(obj, regType)
%CHARTBYTYPE - Return the slice of Charts in this Atlas of the specified RegType
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 20-May-2021; Last revision: 20-May-2021

chartSlice = obj.Chart([obj.Chart.RegType] == regType);
end % end chartbytype

% Revision History:
%{

%}
