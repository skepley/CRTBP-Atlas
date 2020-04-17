classdef Connection < Handle
%CONNECTION - One line description of the class (H1 line)
%
%   CONNECTION() - A more detailed description of the class
%
%   CONNECTION constructor syntax:
%       ConnectionObj = $NAME()
%
%   CONNECTION properties:
%       Property 1 - description
%       Property 2 - description
%
%   CONNECTION methods:
%       Method 1 - description
%       Method 2 - description
%
%   Examples: 
%       Line 1 of example
%       Line 2 of example
%       Line 3 of example
%
%   Subclasses: none
%   Superclasses: none
%   Other classes required: none
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_CLASS_NAME1,  OTHER_CLASS_NAME2

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 18-Jun-2019; Last revision: 18-Jun-2019
% 
%   ToDo: 
%   item1 - 
%   item2 - 


%% ---------------------------------------- Properties ----------------------------------------
    properties
        Chart % A column vector of two Chart objects
        Atlas % A column vector pointing to the atlases containing the Charts
        Local % Local coordinates for the point of intersection
        Global % Global coordinates for the point of intersection
        Phase % A point in phase space where the charts intersect
        Orbit % A k-by-n array whose rows are points in R^n along the connection
        FlightTime % Integration time between both local manifolds
    end % end properties

    
%% ---------------------------------------- Methods ----------------------------------------
    methods
        function obj = Connection(chart1, chart2, atlas1, atlas2, localIntersection, varargin)
        %CONNECTION - class constructor

        end % end class constructor
    end % end methods

end % end classdef

% Revision History:
%{

