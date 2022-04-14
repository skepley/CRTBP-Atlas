%L4_ATLAS - Compute an atlas for the (un)stable manifolds of L4 in the standard CRTBP
%
%   Description:
%       L4_ATLAS description
%
%   Output:
%       L4_ATLAS output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 04-Apr-2019; Last revision: 28-Dec-2019

whichPC = 'mac';
switch whichPC
    case 'mac' % macbook
        computerPath = '/Users/sk2011/';
    case 'lenovo' % ubuntu
        computerPath = '/home/shane/';
end
addpath(genpath([computerPath,'Dropbox/Matlab/IMP']))
addpath([computerPath, 'Dropbox/TaylorDocs/CRFBP_regularization/codes'])
clear all
clc

%% lift local manifolds
isValid = true; % interval or float coefficients
tunstable = 3.0; %...unstable manifold
truncation = [30,10]; % [space, time]
nUnstable = 1;  % number of linear pieces to subdivide local manifold boundary


% lift unstable manifold
load manifoldData_L4_N30
Pu = cat(3,Au,Bu,Cu,Du,Ru,Su);
if ~exist('initUnstable') || ~isequal(nUnstable,length(initUnstable))
    theta = linspace(0,2*pi,nUnstable);
    theta(end+1) = theta(1);
    nodes = exp(1i*theta);
    initUnstable = cell(1,nUnstable);
    initUnstableError = zeros(1,nUnstable);
    % interval enclosure of endpoints not necessary so we enclose after computing exp.
    for k = 1:nUnstable
        parmNode = nodes(k:k+1); % get exp(i*theta) as float and then enclose it.
        parmLine = intval([.5*(parmNode(2) + parmNode(1)), .5*(parmNode(2) - parmNode(1))]); % convert to intval if necessary.
        parmLine(2,:) = conj(parmLine(1,:));
        clear liftUnstable
        for j = 1:6
            unStableError = zeros(1,7);
            thisCord = squeeze(Pu(j,:,:));
            liftCord = real(polycompose(thisCord,parmLine(1,:),parmLine(2,:))); % lift coordinate
            thisTrunc = min([modes(2),length(liftCord)]); % set truncation for shrinkwrapping
            [liftUnstable(j,:),unStableError(j)] = shrinkwrap(liftCord,thisTrunc,rs_validated); % shrinkwrap local and lifting errors into the tail. Improves the numerical stability while integrating.
        end
        initUnstableError(k) = max(unStableError);
        initUnstable{k} = liftUnstable;
    end
end

 



