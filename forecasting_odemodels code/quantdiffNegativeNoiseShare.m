function [share, info] = quantdiffNegativeNoiseShare(fitcurve, scale, distBoot, nSeries)
%QUANTDIFFNEGATIVENOISESHARE  Expected share of bootstrap draws that fall below zero.
%
%   [share, info] = quantdiffNegativeNoiseShare(fitcurve, scale, distBoot, nSeries)
%
%   fitcurve  fitted curve the bootstrap noise is added to (stacked, one block
%             of equal length per fitted series, as in the runners)
%   scale     noise scale used by AddErrorStructure: the standard deviation for
%             normal noise (distBoot = 0), the Laplace scale b for distBoot = 6
%   distBoot  noise model of the bootstrap (dist1_boot in the runners)
%   nSeries   number of fitted series, length(vars.fit_index) (default 1)
%
%   share     mean over the noisy points of P(draw < 0):
%               normal   P = Phi(-mu/sd)        = 0.5*erfc(mu/(sd*sqrt(2)))
%               Laplace  P = 0.5*exp(-mu/b)
%             with mu = |fitted value|, as AddErrorStructure uses. NaN for the
%             other noise models (Poisson and negative binomial draws cannot
%             be negative) and when SCALE is not positive.
%   info      struct: nPoints, minMean, worstPoint (largest single P)
%
%   The first point of every series is left out: the runners reset it to the
%   observed value in each bootstrap dataset, so it carries no noise.
%
%   Why it matters: normal and Laplace noise use ONE scale for the whole
%   series. Where the fitted counts are small compared with that scale (the
%   start and the tail of an epidemic) a large part of the draws is negative.
%   The draws are not truncated, because clipping them at zero would bias the
%   refits upward; the runners warn instead when the share is large, since it
%   means the constant-scale noise model fits the data poorly at low counts.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 4 || isempty(nSeries), nSeries = 1; end

share = NaN;
info = struct('nPoints', 0, 'minMean', NaN, 'worstPoint', NaN);
if ~(isequal(distBoot, 0) || isequal(distBoot, 6))
    return
end
if ~isnumeric(scale) || ~isscalar(scale) || ~isfinite(scale) || scale <= 0
    return
end

mu = abs(double(fitcurve(:)));
if mod(numel(mu), nSeries) ~= 0
    error('quantdiffNegativeNoiseShare:series', ...
        '%d fitted values cannot be split into %d series of equal length.', numel(mu), nSeries);
end
nT = numel(mu) / nSeries;
noisy = true(size(mu));
noisy(1:nT:end) = false;                 % first point of every series: reset to the data, no noise
noisy = noisy & isfinite(mu);
if ~any(noisy)
    return
end
mu = mu(noisy);

if distBoot == 0
    p = 0.5 * erfc(mu ./ (scale * sqrt(2)));     % normal: Phi(-mu/sd)
else
    p = 0.5 * exp(-mu ./ scale);                 % Laplace
end

share = mean(p);
info = struct('nPoints', numel(mu), 'minMean', min(mu), 'worstPoint', max(p));
end
