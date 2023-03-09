%CHECK_REGCRTBP - testing, validation, and debugging for RegCRTBPChart class for 0 dimensional initial data

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 22-Apr-2020; Last revision: 22-Apr-2020

warning('This has not been updated to the new calling syntax for RegCRTBPCharts which requires an energy')
clearvars
close all
clc
computerPath = pcpath('mac');
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath(genpath([computerPath, 'Dropbox/Regularisation3bp/CRTBP Atlas']))
%% ================================================== MAKE SOME CHOICES ==================================================
% Integrator parameters
M = 25; % temporal truncation
truncation = M; % set truncation as a single float
basis = 'Taylor';
mu = .25; % small mass primary
isValid = false;
bdCheck = @(obj, boundaryChart, maxTau)boundarycheck(obj, boundaryChart, maxTau);
odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);

% Integrator time stepping
initialTime = 0;
tauGuess = 0.1; % initial timestep to scale each Taylor step
maxTau = .3;
tf = linspace(0, maxTau, 250); % plotting time points
regType = 0;

switch regType
    case 0  % =============== CHECK f_0 AGAINST RUNGE KUTTA ===============
        parameter = mu;
        % set up some initial data at a point
        p1 = [0.7104, 0.2973, 0.2554, 0.2068];
        initialData = p1.';
        % Integrate multiple timesteps
%         bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
        cht = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess);
%         % Integrate multiple timesteps
%         A0 = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
%         while ~isempty(A0.LeafStack)
%             A0.growboundary()
%         end
        
        % plot orbits
        figure;
        hold on      
        try
            ob = cht.eval(tf');
            ob = cell2mat(ob')';
        catch
            
        end
        
        plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',2)
        
        % set up rk45 evaluation
        VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
        initCondition = ob(1,:)';
        [~,sol] = ode45(VF,tf,initCondition,odeOptions);
        plot3(sol(:,1),sol(:,3),sol(:,2),'go')
        view(-50,50)
        
    case 1  % =============== CHECK f_1 AGAINST RUNGE KUTTA ===============
        C = 3;
        parameter = [mu,C];
        
        % set up some initial data on a line
        p1 = [0.7104, 0.2973, 0.2554, 0.2068];
        initialData = p1.';
        bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
        cht = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess);
        
%         % Integrate multiple timesteps
%         bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
%         maxTau = .3;
%         A1 = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
%         while ~isempty(A1.LeafStack)
%             A1.growboundary()
%         end
        
        % plot orbits
        figure;
        hold on
        s = 0;  

        ob = cht.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']);
        ob = cell2mat(ob')';
        plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',2)
        
        % set up rk45 evaluation
        VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
        initCondition = ob(1,:)';
        [~,sol] = ode45(VF,tf,initCondition,odeOptions);
        plot3(sol(:,1),sol(:,3),sol(:,2),'go')
        view(-50,50)
        
    case 2  % =============== CHECK f_2 AGAINST RUNGE KUTTA ===============
        C = 3;
        parameter = [mu,C];
        % set up some initial data on a circle
        thetaMid = 0;
        thetaRadius = .0001;
        expCoeff = sqrt(8*mu)*exp(1i*thetaMid)*((1i*thetaRadius).^(0:N-1)./factorial(0:N-1));
        xInitial = zeros(1,N);
        pInitial = real(expCoeff);
        yInitial = zeros(1,N);
        qInitial = imag(expCoeff);
        initialData = [xInitial;pInitial;yInitial;qInitial];
        
        % Integrate multiple timesteps
        bd = RegCRTBPChart(initialData, basis, initialTime, truncation, parameter, regType, 'InitialScaling', tauGuess, 'boundary', true);
        maxTau = .3;
        A2 = RegCRTBPAtlas(bd, tauGuess, bdCheck, @advectioncheck, 'MaxTau', maxTau);
        while ~isempty(A2.LeafStack)
            A2.growboundary()
        end
        
        % plot orbits
        figure;
        hold on
        s = 0;
        ob = mid(A2.orbit(s, tf));
        plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',2)
        
        % set up rk45 evaluation
        VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
        initCondition = ob(1,:)';
        [~,sol] = ode45(VF,tf,initCondition,odeOptions);
        plot3(sol(:,1),sol(:,3),sol(:,2),'go')
        view(-50,50)
end % switch

%% ================================================== PLOT AGAINST RK45 ==================================================
%         % close all
%         % figure;
%         hold on
%
%         % plot orbits
%         tf = linspace(0, maxTau, 250);
%         s = 0;
%         try
%             ob = cell2mat(orbit.eval([s*ones(length(tf),1),linspace(0,1,length(tf))']));
%         catch
%             ob = mid(A.orbit(s, tf));
%         end
%         plot3(ob(:,1), ob(:,3), ob(:,2), 'k', 'LineWidth',1)
%
%         % d1 = sum(abs(ob(:,5) - 1./sqrt((ob(:,1) - mu).^2 + ob(:,3).^2)))
%         % d2 = sum(abs(ob(:,6) - 1./sqrt((ob(:,1)+1-mu).^2 + ob(:,3).^2)))
%
%         % set up rk45 evaluation
%         odeOptions = odeset('RelTol',1e-13,'AbsTol',1e-13);
%         VF = @(t,x)rk45regvectorfield(t,x,parameter,regType);
%         initCondition = ob(1,:)';
%         [~,sol] = ode45(VF,tf,initCondition,odeOptions);
%         plot3(sol(:,1),sol(:,3),sol(:,2),'go')
%         view(-50,50)
%         dealfig()
%