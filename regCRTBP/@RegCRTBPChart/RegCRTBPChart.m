classdef RegCRTBPChart < Chart
    %REGCRTBPCHART - A parameterization of a single Taylor timestep for the regularized circular restricted three body problem (CRTBP).
    %
    %   REGCRTBPCHART constructor syntax:
    %       RegCRTBPChartObj = $NAME(initialData, basis, initialTime, truncation, parameter, regType, varargin) constructs an instance of this class.
    %           initialData: A k-by-d, array of initial data. For k = 4 the automatic differential is automatically computed and for k = 5,6 the extra
    %               dimensions specify the automatic differentiation coordinates for F1/F2 and F0 fields respectively.
    %           basis - One of the following: 'Taylor', {'Taylor','Taylor'}, {'Chebyshev','Taylor'}. The first 2 are equivalent.
    %           initialTime - For keeping track of a global time variable. It has no role in the integration since all 3 fields are autonomous.
    %           parameter - Either mu or [mu, C] for standard and regularized fields respectively.
    %           regType - A value in {0,1,2} for standard (F0) or regularized (F1/F2) fields.
    %           varags:
    %
    %   REGCRTBPCHART properties:
    %       RegType - A value in {0,1,2} which defines which coordinates this chart is computed with respect to. F0 for standard CRTBP vector field,
    %           F1 for the regularized vector field desingularized at the large primary, and F2 for the regularized field desingularized at the small primary.
    %       Parameter - Specify the mass of the small primary by mu. For the regularized field, this has a second value which specifies the energy of regularization.
    %
    %   REGCRTBPCHART methods:
    %       All methods have been moved out of the class definition file. Their documentation can be found in their respective files.
    %
    %   Subclasses: none
    %   Superclasses: Chart
    %   Other classes required: Scalar
    %   Other m-files required: IntLab
    %   MAT-files required: none
    
    %   Author: Shane Kepley
    %   email: shane.kepley@rutgers.edu
    %   Date: 07-Mar-2019; Last revision: 20-Feb-2020
    %
    %   ToDo: 1. RegCRTBPChart/hotswap needs to inherit the initial error when coordinate swapping before this computation can be made rigorous. 
    
    properties
        RegType; % 0, 1, or 2 denotes a chart with respect to F0, F1, or F2.
        Parameter; % mu or [mu,C]
        LocalTime; % Keep track of the regularized local time so this information can be passed to the BVP solver.
    end
    
    properties(Hidden = true)
    end
    
    %% ---------------------------------------METHODS--------------------------------------------------%%
    methods
        function obj = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, varargin)
            % class constructor
            if(nargin > 0)
                % parse input
                p = inputParser;
                addRequired(p,'initialData') % 4 or 5-by-spatialDegree array of double or intval
                addRequired(p,'basis') % Taylor, Fourier, or Chebyshev
                addRequired(p,'initialTime') % 1-by-1 double
                addRequired(p,'truncation') % 1-by-(1+d) integer
                addRequired(p,'parameter') % [mu,C] mu is the small mass, C is the regularization energy
                addRequired(p,'regType') % regularization type of this chart
                addParameter(p,'MaxTau', Inf) % scalar double
                addParameter(p,'FixTau', false) % logical
                addParameter(p,'InitialScaling', 1) % positive scalar double
                addParameter(p,'InitialError', 0) % % ell^1 error bounds on initial data.
                addParameter(p,'SpatialSpan', [-1,1]) % 1-by-2 double
                addParameter(p,'ParentId',[]) % RegCRTBPChart from the previous timestep
                addParameter(p,'Boundary', false)
                % addParameter(p,'Weight', 'ones') % vector of weights for the ell_1 space to work in
                
                % parse variable arguments
                parse(p, initialData, basis, initialTime, truncation, parameter, regType, varargin{:})
                obj.RegType = p.Results.regType; % set the regularization type
                obj.TimeSpan = p.Results.initialTime;
                obj.Tau = p.Results.InitialScaling; % Setting << 1 improves precision of multiplication when FFT is required. The recursion solves instead x' = Lf(x,Lt).
                obj.InitialError = p.Results.InitialError;
                obj.SpatialSpan = p.Results.SpatialSpan;
                obj.ParentHandle = p.Results.ParentId;
                % obj.Weight = p.Results.Weight;
                isBoundary = p.Results.Boundary; % generate the initial data as a boundary chart but don't advect it.
                
                % set properties for instance of RegCRTBPChart
                obj.Coordinate = Scalar; % initialize coordinates
                obj.NumericalClass = class(initialData); % intval or double coefficients
                obj.Truncation = truncation;
                
                % phase space dimension
                if isequal(obj.RegType, 0)
                    phaseDimension = 6;
                else
                    phaseDimension = 5;
                end
                obj.Dimension = [length(truncation)-1, phaseDimension]; % [chart dimension, phase dimension]
                obj.ClassConstructor = @RegCRTBPChart; % handle to create more charts of this subclass
                
                % set parameters to intervals when doing validated computations
                if isequal(obj.NumericalClass,'intval')
                    obj.Parameter = intval(p.Results.parameter);
                    obj.IsValid = true;
                elseif isequal(obj.NumericalClass,'double')
                    obj.Parameter = p.Results.parameter;
                    obj.IsValid = false;
                end
                
                % automatic differentiation of coordinates
                obj.setinitialdata(initialData, basis);  % write initial data to Scalar
                if size(initialData,1) == 4 % need to compute validated initial data for r
                    obj.validateinitialdata(); % validated computation of automatic differentiation coordinates
                end
                
                % initialize coordinates as Scalars
                if isequal(obj.Dimension(1), 0)
                    obj.Coordinate = obj.InitialData.deepcopy();
                    % initial data has a dimension = 0 and truncation = 1 as required
                else
                    obj.Coordinate = Scalar(obj.InitialData, basis, obj.Truncation(2:end));
                end
                obj.Basis = obj.Coordinate(1).Basis; % set basis for chart as a cell array
                
                % advect this initial data
                if ~isBoundary
                    obj.advect(obj.Tau); % compute chart coefficients for advected image
                end
            end % end if nargin
        end % end class constructor
    end % end methods
end % end RegCRTBPChart


% Revision History:
%{
06-Apr-2019 - Merged the classes for all 3 regularization types into a single class. The CRTBPChart class is deprecated.
23-May-2019 - Added support for boundary charts in Chebyshev basis with integration in Taylor basis. Tested on local (un)stable manifolds for Lyapunov
    orbits for L1,L2,L3 provided by Maxime Murray.
%}
