function [nEff, info] = quantdiffEffectiveN(nObs, nSeries, fixI0)
%QUANTDIFFEFFECTIVEN  Number of observations that the fit actually has to explain.
%
%   [nEff, info] = quantdiffEffectiveN(nObs, nSeries, fixI0)
%
%   nObs     observations in the calibration window, all fitted series stacked
%            (length(ydata) = window length x number of fitted series)
%   nSeries  number of fitted series, length(vars.fit_index)
%   fixI0    params.fixI0
%
%   With params.fixI0 = 1 the initial value of every fitted state is set to the
%   first observation of its series, and the first point of the fitted curve IS
%   that initial value (quantdiffObservationCurve), so the first residual of
%   each series is zero by construction: that point carries no information
%   about the parameters. It is therefore not counted,
%
%       nEff = nObs - nSeries          (fixI0 = 1)
%       nEff = nObs                    (fixI0 = 0)
%
%   and nEff replaces n wherever the sample size enters: AICc (the n*log(SSE/n)
%   and 2n*log(2*SAD/n) terms and the small-sample correction), the residual
%   variance of the normal bootstrap noise, the profile deviance of methods 0
%   and 6, and the "enough data points" checks. The objective value itself is
%   not changed (a zero residual adds nothing to a sum of squares or of
%   absolute deviations).
%
%   With fixI0 = 0 nothing is subtracted: the free initial values are counted
%   as estimated quantities (get_nparams), so n - k is already right, and
%   whether an initial value only matches the first point or also drives the
%   dynamics depends on the model, which the options file does not say.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if ~isnumeric(nObs) || ~isscalar(nObs) || ~isfinite(nObs) || nObs < 1 || nObs ~= fix(nObs)
    error('quantdiffEffectiveN:nObs', 'nObs must be a positive integer.');
end
if ~isnumeric(nSeries) || ~isscalar(nSeries) || ~isfinite(nSeries) || nSeries < 1 || nSeries ~= fix(nSeries)
    error('quantdiffEffectiveN:nSeries', 'nSeries must be a positive integer.');
end
if mod(nObs, nSeries) ~= 0
    error('quantdiffEffectiveN:series', ...
        '%d observations cannot be split into %d series of equal length.', nObs, nSeries);
end

pinned = 0;
if ~isempty(fixI0) && fixI0 == 1
    pinned = nSeries;
end
nEff = nObs - pinned;
info = struct('nObs', nObs, 'nSeries', nSeries, 'pinnedPoints', pinned);
end
