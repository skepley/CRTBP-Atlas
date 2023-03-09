function isIntersecting = intersectellonebox(chart1, chart2)
%INTERSECTELLONEBOX - Checks if these two charts intersect.
%
%   Syntax:
%       isIntersecting = INTERSECTELLONEBOX(chart1, chart2) is false if this pair of charts is rigorously nonintersecting via ell-one boxes computed
%       using IntLab. Otherwise returns true implying Newton's method is required.
%
%   Inputs:
%       chart1 - A RegCRTBPChart
%       chart2 - A RegCRTBPChart
%
%   Outputs:
%       isIntersecting - true or false
%
%   Subfunctions: none
%   Classes required: RegCRTBPChart
%   Other m-files required: IntLab
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 12-Jun-2019; Last revision: 17-May-2021

% first rule out using ell-1 box method
box1 = chart1.ellonebox();
box2 = chart2.ellonebox();
mu = chart1.Parameter(1);

if ~isequal(chart1.RegType, chart2.RegType) % charts are not in same coordinates
    if ~isequal(chart1.RegType, 0) % chart1 is not in F0 coordinates
        box1 = CRTBP2reg(box1(1:4), mu, -1*chart1.RegType); % change chart 1 from F1/F2 into F0 coordinates
    end
    
    if ~isequal(chart2.RegType, 0) % chart2 is not in F0 coordinates
        box2 = CRTBP2reg(box2(1:4), mu, -1*chart2.RegType); % change chart 2 from F1/F2 into F0 coordinates
    end
    
end

ellOneDistance = abs(box1(1:4) - box2(1:4));
isIntersecting = ~any(inf(ellOneDistance) > 0); % returns false if connection ruled out by ell^1 box
end % end intersectellonebox

% Revision History:
%{
17 May 2021 - Changed intersection of charts in different RegType to apply CRTBP2Reg function to the ellone box to
obtain interval enclosures of images in F0 coordinates.
%}
