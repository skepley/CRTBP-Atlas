classdef RegCRTBPConnection < handle
    %REGCRTBPCONNECTION - One line description of the class (H1 line)
    %
    %   REGCRTBPCONNECTION() - Connection class for the GCS algorithm for the CRTBP with Levi-Cevita regularization.
    %
    %   REGCRTBPCONNECTION constructor syntax:
    %       RegCRTBPConnectionObj = $NAME()
    %
    %   REGCRTBPCONNECTION properties:
    %       Property 1 - description
    %       Property 2 - description
    %
    %   REGCRTBPCONNECTION methods:
    %       Method 1 - description
    %       Method 2 - description
    %
    %   Examples:
    %       Line 1 of example
    %       Line 2 of example
    %       Line 3 of example
    %
    %   Subfunctions: none
    %   Classes required: none
    %   Other m-files required: none
    %   MAT-files required: none
    
    %   Author: Shane Kepley
    %   email: s.kepley@vu.nl
    %   Date: 23-Feb-2023;
    %
    %   ToDo:
    %   item1 -
    %   item2 -
    
    
    %% ---------------------------------------- Properties ----------------------------------------
    properties
        StableChart
        UnstableChart
        ConnectionEnergy
        LocalIntersection  % [stableSpace, stableTime, unstableSpace] coordinates for the interesection in the local coordinate system which is a subset of D^d
        GlobalIntersection  % [stableSpace, stableTime, unstableSpace, unstableTime] = global coordinates for the intersection:
        ConnectionTime  % time of flight in physical time between source manifold and target manifold
        LocalMaps  % [stableLocalMap, unstableLocalMap]
        TrueOrbit  % heuristics for deciding if the orbit is a connection or pseudo connection
        Parameter % mass ratio of the small primary
        
        
    end % end properties
    
    
    %% ---------------------------------------- Methods ----------------------------------------
    methods
        function obj = RegCRTBPConnection(connectionIntersectionData, stableLocalMap, unstableLocalMap)
            %REGCRTBPCONNECTION - class constructor
            if nargin > 0
                warning('OFF', 'Chart:eval');  % turn off Chart.eval warning messages
                warning('OFF', 'Chart:local2global');  % turn off Chart.eval warning messages
                
                obj.StableChart = connectionIntersectionData{1};
                obj.UnstableChart = connectionIntersectionData{2};
                obj.ConnectionEnergy = obj.StableChart.RegEnergy;
                obj.LocalIntersection = [connectionIntersectionData{3}; 0];  % [stableSpace; stableRegTime; unstableSpace; unstableTime]
                
                localStablePhysTime = taylorregtime(obj.StableChart, obj.LocalIntersection(2));
                
                obj.GlobalIntersection = [obj.StableChart.local2global([obj.LocalIntersection(1), localStablePhysTime]),...
                    obj.UnstableChart.local2global([obj.LocalIntersection(3), 0])].'; % [stableSpace, stableTime, unstableSpace, unstableTime]
                obj.ConnectionTime = obj.GlobalIntersection(4) - obj.GlobalIntersection(2);  % time of flight in physical time
                obj.TrueOrbit = connectionIntersectionData{4};
                obj.LocalMaps = {stableLocalMap, unstableLocalMap};
                obj.Parameter = obj.StableChart.Parameter;
            end
        end % end class constructor
    end % end methods
end
