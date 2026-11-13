function [b, info] = estimateLaplaceScale(ydata, fitcurve, numModelParams, nSeries, fixI0)
%ESTIMATELAPLACESCALE  Laplace scale of the data around the fitted curve (SAD fits).
%
%   [b, info] = estimateLaplaceScale(ydata, fitcurve, numModelParams, nSeries, fixI0)
%
%   Scale of the Laplace bootstrap noise (AddErrorStructure, dist1 = 6) for
%   models fitted by the sum of absolute deviations (method1 = 6):
%
%       b = sum( |y_i - mu_i| ) / (n_used - numModelParams/2)
%
%   i.e. the maximum-likelihood Laplace scale SAD/n (the value behind getAICc,
%   case 6, and the profile likelihood) with a small-sample correction. For an
%   L1 fit with p parameters the expected reduction of the sum of absolute
%   residuals is p/(4 f(0)), with f(0) the error density at zero; for Laplace
%   errors f(0) = 1/(2b), so E[SAD] ~ b (n - p/2). Subtracting the full
%   parameter count, as for least squares, over-corrects. Simulations of SAD
%   fits of a 3-parameter growth model agree: n - p/2 gave a bias of +1% to +2%
%   for n = 19 and 39, n - p gave +5% to +12%, and plain n gave -3% to -6%.
%   Unlike the former SD/sqrt(2), the estimate does not square the residuals,
%   so a few outliers inflate it far less.
%
%   ydata, fitcurve  stacked observations and fitted curve (same length; one
%                    block of equal length per fitted series, as in the runners)
%   numModelParams   estimated quantities of the curve (model parameters and
%                    free initial values), NOT counting the Laplace scale itself;
%                    the runners pass numparams - 1
%   nSeries          number of fitted series stacked in ydata (default 1)
%   fixI0            1 -> the first point of each series is fixed to the data
%                    by construction (zero residual) and is excluded (default 0)
%
%   b     the estimate; NaN when no more than numModelParams points are usable
%         (the caller decides what to do, as with estimatePearsonDispersion)
%   info  struct: nUsed, dof, excludedFixed
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 4 || isempty(nSeries), nSeries = 1; end
if nargin < 5 || isempty(fixI0),   fixI0 = 0;   end

y  = double(ydata(:));
mu = double(fitcurve(:));
if numel(y) ~= numel(mu)
    error('estimateLaplaceScale:size', ...
        'ydata (%d) and fitcurve (%d) must have the same number of elements.', numel(y), numel(mu));
end
if mod(numel(y), nSeries) ~= 0
    error('estimateLaplaceScale:series', ...
        '%d observations cannot be split into %d series of equal length.', numel(y), nSeries);
end
nT = numel(y) / nSeries;

use = isfinite(y) & isfinite(mu);
fixedMask = false(size(y));
if fixI0
    fixedMask(1:nT:end) = true;          % first point of every series
end
use = use & ~fixedMask;

dof = nnz(use) - numModelParams/2;
if nnz(use) <= numModelParams
    b = NaN;
else
    b = sum(abs(y(use) - mu(use))) / dof;
end

info = struct('nUsed', nnz(use), 'dof', dof, 'excludedFixed', nnz(fixedMask));
end
