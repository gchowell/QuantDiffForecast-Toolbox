function [perfC, perfF] = computeScoresFromQuantiles(Q, levels, datalatest, nCalib, horizon)
%COMPUTESCORESFROMQUANTILES  MAE, MSE, 95% coverage and WIS of a forecast given by its quantiles.
%
%   [perfC, perfF] = computeScoresFromQuantiles(Q, levels, datalatest, nCalib, horizon)
%
%   Q           (nCalib + horizon) x numel(levels) matrix of predictive
%               quantiles, one row per time point (calibration rows first)
%   levels      probability levels of the columns of Q. They must include
%               the 23 levels used by computeWIS:
%               [0.01 0.025 0.05:0.05:0.95 0.975 0.99]
%   datalatest  observed series, [time value], starting at the first
%               calibration row (as saved by the runners)
%   nCalib      number of calibration rows
%   horizon     number of forecast steps
%
%   perfC   1 x 4   [MAE MSE Coverage95 WIS] over the calibration rows
%   perfF   horizon x 5   [h MAE MSE Coverage95 WIS]; row h averages forecast
%           steps 1..h, the convention of computeforecastperformance and
%           computeWIS. Empty when horizon = 0 or the observations do not
%           cover the full horizon (the same rule as those functions).
%
%   The scores follow computeforecastperformance and computeWIS exactly, except
%   that the quantiles are supplied instead of being taken from draws:
%     point forecast  the 0.5 quantile
%     coverage        share of observations inside [q0.025, q0.975], in percent
%     WIS             (1/(K + 1/2)) * ( 0.5*|y - median| + sum_k (alpha_k/2) * IS_k ),
%                     alpha = [0.02 0.05 0.1:0.1:0.9], K = 11
%   Interval bounds below zero are set to zero, as in those functions.
%
%   Purpose: a linear-pool ensemble has EXACT quantiles (the weighted quantiles
%   of all members' draws, getensemble_linearpool_from_curves). Scoring those
%   removes the sampling noise of scoring a random sample of trajectories drawn
%   from the pool (about 1% of the WIS with 6000 draws per member, about 5%
%   with 400). run_ensemble scores the ensemble, its members and the baselines
%   through this function, so they are compared on the same footing.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

alphas = [0.02 0.05 0.1:0.1:0.9];
K = numel(alphas);
levels = double(levels(:)');

if size(Q, 2) ~= numel(levels)
    error('computeScoresFromQuantiles:size', ...
        'Q has %d columns but %d levels were given.', size(Q, 2), numel(levels));
end
if size(Q, 1) < nCalib + horizon
    error('computeScoresFromQuantiles:rows', ...
        'Q has %d rows; %d calibration rows + %d forecast steps are needed.', size(Q, 1), nCalib, horizon);
end

% Columns of the levels the scores need (matched with a tolerance, since
% 0.05:0.05:0.95 is not exact in floating point).
iMed = localColumn(levels, 0.5);
iLo  = zeros(1, K);  iHi = zeros(1, K);
for k = 1:K
    iLo(k) = localColumn(levels, alphas(k) / 2);
    iHi(k) = localColumn(levels, 1 - alphas(k) / 2);
end
iL95 = localColumn(levels, 0.025);
iU95 = localColumn(levels, 0.975);

y = datalatest(:, 2);
nObs = numel(y);
if nObs < nCalib
    error('computeScoresFromQuantiles:data', ...
        'The observed series has %d rows; %d calibration rows are needed.', nObs, nCalib);
end

% Per-row scores for every row that has an observation.
nRows = min(nObs, nCalib + horizon);
absErr = nan(nRows, 1);  sqErr = nan(nRows, 1);  inside = nan(nRows, 1);  wis = nan(nRows, 1);
for t = 1:nRows
    m = Q(t, iMed);
    absErr(t) = abs(y(t) - m);
    sqErr(t)  = (y(t) - m)^2;
    L95 = max(Q(t, iL95), 0);  U95 = max(Q(t, iU95), 0);
    inside(t) = double(y(t) >= L95 && y(t) <= U95);
    s = 0;
    for k = 1:K
        a = alphas(k);
        L = max(Q(t, iLo(k)), 0);  U = max(Q(t, iHi(k)), 0);
        IS = (U - L) + (2 / a) * (L - y(t)) * (y(t) < L) + (2 / a) * (y(t) - U) * (y(t) > U);
        s = s + (a / 2) * IS;
    end
    wis(t) = (0.5 * absErr(t) + s) / (K + 0.5);
end

c = 1:nCalib;
perfC = [mean(absErr(c)) mean(sqErr(c)) 100 * mean(inside(c)) mean(wis(c))];

perfF = [];
if horizon > 0 && nObs >= nCalib + horizon
    perfF = zeros(horizon, 5);
    for h = 1:horizon
        f = nCalib + (1:h);
        perfF(h, :) = [h mean(absErr(f)) mean(sqErr(f)) 100 * mean(inside(f)) mean(wis(f))];
    end
end
end


function j = localColumn(levels, p)
% Column of LEVELS equal to P (tolerance 1e-9); error if it is missing.
j = find(abs(levels - p) < 1e-9, 1, 'first');
if isempty(j)
    error('computeScoresFromQuantiles:level', ...
        'The quantile level %g is required but is not among the supplied levels.', p);
end
end
