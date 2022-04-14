classdef CRTBPChart < Chart
    %CRTBPCHART - Validated Taylor timestep for the standard circular restricted three body problem.
    %
    %   CRTBPCHART constructor syntax:
    %       CRTBPChartObj = $NAME()
    %
    %   CRTBPCHART properties:
    %       Property 1 - description
    %       Property 2 - description
    %
    %   CRTBPCHART methods:
    %       Method 1 - description
    %       Method 2 - description
    %
    %   Subclasses: none
    %   Superclasses: Chart
    %   Other classes required: Scalar, Chart
    %   Other m-files required: none
    %   MAT-files required: none
    
    %   Author: Shane Kepley
    %   email: shane.kepley@rutgers.edu
    %   Date: 20-Mar-2019; Last revision: 20-Mar-2019
    %
    %   ToDo:
    %   item1 -
    %   item2 -
    
    
    properties
        Parameter; % [mu]
    end
    
    properties(Hidden = true)
    end
    
    %% ---------------------------------------METHODS--------------------------------------------------%%
    methods
        function obj = CRTBPChart(initialData, basis, initialTime, truncation, parameter, varargin)
            % class constructor
            obj.ClassConstructor = @CRTBPChart; % handle to create more charts of this subclass
            if(nargin > 0)
                % parse input
                p = inputParser;
                addRequired(p,'initialData') % % 4 or 6-by-spatialDegree array of double or intval
                addRequired(p,'basis') % Taylor, Fourier, or Chebyshev
                addRequired(p,'initialTime') % 1-by-1 double
                addRequired(p,'truncation') % 1-by-(1+d) integer
                addRequired(p,'parameter') % [mu1, mu2]
                addParameter(p,'Direction', 1) % 1 or -1 (forward or backward integration)
                addParameter(p,'MaxTau', Inf) % scalar double
                addParameter(p,'FixTau', false) % logical
                addParameter(p,'InitialScaling', 1) % positive scalar double
                addParameter(p,'InitialError', 0) % % ell^1 error bounds on initial data.
                addParameter(p,'SpatialSpan', [-1,1]) % 1-by-2 double
                addParameter(p,'ParentId',[]) % CRTBPChart from the previous timestep
                addParameter(p,'Boundary', false)
                
                % parse variable arguments
                parse(p, initialData, basis, initialTime, truncation, parameter, varargin{:})
                obj.TimeSpan = p.Results.initialTime;
                obj.Tau = p.Results.InitialScaling; % Setting << 1 improves precision of multiplication when FFT is required. The recursion solves instead x' = Lf(x,Lt).
                obj.InitialError = p.Results.InitialError;
                obj.SpatialSpan = p.Results.SpatialSpan;
                obj.ParentHandle = p.Results.ParentId;
                isBoundary = p.Results.Boundary; % generate the initial data as a boundary chart but don't advect it.
                
                % set properties for instance of CRTBPChart
                obj.Coordinate = Scalar; % initialize coordinates
                obj.NumericalClass = class(initialData); % intval or double coefficients
                obj.Truncation = truncation;
                obj.Dimension = [length(truncation)-1, 6]; % [chart dimension, phase dimension]
                
                % automatic differentiation of 5th and 6th coordinates
                if size(initialData,1) == 4 % need to compute validated initial data for r
                    obj.setinitialdata(initialData, basis);  % write initial data for first 4 coordinates to Scalar
                    obj.validateinitialdata(obj);
                else % all 6 coordinates of initial data are provided
                    obj.setinitialdata(initialData, basis);  % write initial data for all 6 coordinates to Scalar
                end
                
                % set mass parameters
                if isequal(obj.NumericalClass,'intval')
                    obj.Parameter = intval(p.Results.parameter);
                    obj.IsValid = true;
                elseif isequal(obj.NumericalClass,'double')
                    obj.Parameter = p.Results.parameter;
                    obj.IsValid = false;
                end
                
                % initialize coordinates as Scalars
                if isequal(obj.Dimension(1), 0)
                    obj.Coordinate = obj.InitialData.deepcopy();
                    % initial data has a dimension = 0 and truncation = 1 as required
                else
                    obj.Coordinate = Scalar(obj.InitialData, 'Taylor', obj.Truncation(2:end));
                end
                if ~isBoundary
                    obj.advect(obj.Tau); % compute chart coefficients for advected image
                end
            end % end if nargin
        end % end class constructor
        
        function advect(obj, tau)
            %ADVECT Advect the (d-1) dimensional boundary of a chart to generate the d dimensional interior chart
            
            % Increase Chart dimension and Scalar coordinate dimensions by 1 to prepare for new coefficients to be appended
            obj.Dimension(1) = obj.Dimension(1) + 1; % increase domain dimension for Chart
            for iCoordinate = 1:obj.Dimension(2) % update dimensions of Scalar coordinates
                obj.Coordinate(iCoordinate).Dimension = obj.Dimension(1);
                obj.Coordinate(iCoordinate).Truncation = [1, obj.Truncation(2:end)];
            end
            obj.generatecoefficient(tau) % advect this boundary Chart
            obj.Tau = tau; % set step size
            obj.TimeSpan = [obj.TimeSpan, obj.TimeSpan+tau]; % set time span
        end
    end % end methods
end % end CRTBPChart

% Revision History:

%{

%}
