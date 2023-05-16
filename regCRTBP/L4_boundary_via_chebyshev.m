function bdCharts = L4_boundary_via_chebyshev(sourceFile, truncation, type, energy, tauGuess)
%L4_TO_CHEB - Evaluate the boundary of a local invariant manifold attached to L4 using Chebyshev interpolation
%
%   Syntax:
%
% Inputs:
%   - sourceFile: The file containing local manifold data.
%   - type: Type of the system ('stable' or 'unstable').
%   - truncation: Truncation limits as a vector.
%   - energy: Energy value of L4
%   - mu: mass parameter value of the system.
%   - tauGuess: Initial time scaling value to use when integrating.
%
% Outputs:
%   - bdCharts: An array of regCRTBPCharts which parameterize the boundary of the invariant manifold
%
%   Subfunctions: chebyshev_interpolate
%   Classes required: Scalar, regCRTBPChart
%   Other m-files required: Chebfun
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 26-Mar-2023;

% import local manifold data and make into Scalars
basis = {'Taylor', 'Taylor'};
load(sourceFile);
switch type
    case 'stable'
        P = [Scalar(As, basis); Scalar(Bs, basis); Scalar(Cs, basis); Scalar(Ds, basis); Scalar(Rs, basis); Scalar(Ss, basis)];
    case 'unstable'
        P = [Scalar(Au, basis); Scalar(Bu, basis); Scalar(Cu, basis); Scalar(Du, basis); Scalar(Ru, basis); Scalar(Su, basis)];
end

% Interpolate the full boundary in one chart using chebfun to choose the number of modes automatically
chebBd = chebyshev_interpolate(P);
tmpTrunc = max([chebBd.Truncation]);
for idx = length(chebBd):-1:1
    chebCoefs(idx, :) = Scalar.embed(chebBd(idx).Coefficient, tmpTrunc);
end
fullBdChart = RegCRTBPChart(chebCoefs, 'Taylor', 0, [truncation(1), tmpTrunc], energy, mu, 0, 'InitialScaling', tauGuess, 'boundary', true);

bdCharts = fullBdChart.controltailratio(truncation(2), 1e-6); % remesh until most of the ell^1 weight is in the first N coefficients
for bd_idx = 1:length(bdCharts)  % truncate the coefficients of the remeshed charts to order N
    bd = bdCharts(bd_idx);
    for sc_idx = 1:length(bd.Coordinate)
        bd.Coordinate(sc_idx).Coefficient = Scalar.embed(real(bd.Coordinate(sc_idx).Coefficient), truncation(2));
        bd.Coordinate(sc_idx).Truncation = truncation(2);
    end
    bd.Truncation = truncation;
end
end % end L4_to_cheb

function fout = chebyshev_interpolate(P)
    % chebyshev_interpolate performs Chebyshev interpolation on the given data.
    % Inputs:
    % - P: Coefficients for interpolation.
    %
    % Classes required: Scalar
    N = P(1).Truncation(1) - 1;
    for j = length(P):-1:1
        Pj = @(s)P(j).eval([exp(1i*pi*s), exp(-1i*pi*s)]);
        Fj = chebfun(Pj);
        fout(j) = Scalar(fliplr(poly(Fj)), 'Taylor');
    end
    fout = fout';
end
