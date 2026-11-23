function S = quantdiffForecastScores(draws, datalatest, nCalib, horizon, varargin)
%QUANTDIFFFORECASTSCORES  Scores of a forecast for every horizon on its own.
%
%   S = quantdiffForecastScores(draws, datalatest, nCalib, horizon)
%   S = quantdiffForecastScores(..., 'Aggregate', 7)
%   S = quantdiffForecastScores(..., 'BaselineDraws', baseDraws)
%
%   draws       (nCalib + horizon) x M matrix of predictive draws, laid out like
%               forecast2 in the runners: calibration rows first, one column per
%               trajectory (bootstrap curve plus observation noise)
%   datalatest  observed series, [time value], starting at the first
%               calibration row (as saved by the runners)
%   nCalib      number of calibration rows
%   horizon     number of forecast steps
%
%   'Aggregate'      k (default 1). With k > 1 the target is the TOTAL over
%                    blocks of k consecutive forecast steps (k = 7 turns daily
%                    counts into weekly totals): every trajectory is summed
%                    over each block before scoring, and so are the
%                    observations. A trailing incomplete block is dropped.
%   'BaselineDraws'  draws of a reference forecast in the same layout (for
%                    example from baselineForecast); adds its WIS and the
%                    relative WIS to the output
%
%   Output S (struct)
%     S.perHorizon  table, one row per horizon (or block), each scored on its own:
%                     horizon, step_from, step_to   forecast steps covered
%                     time_from, time_to            their time labels (NaN beyond the data)
%                     observed, median              observation and predictive median
%                     abs_error, sq_error           of the median
%                     cover50, cover95              1 if the observation is inside the
%                                                   50% / 95% prediction interval, else 0
%                     WIS                           weighted interval score
%                     WIS_log                       WIS of log(x + 1): errors are relative,
%                                                   so the largest counts do not dominate
%                     WIS_baseline, relWIS          with 'BaselineDraws': the baseline's WIS
%                                                   and WIS ./ WIS_baseline
%     S.cumulative  table, row h = mean over horizons 1..h of the rows above
%                   (MAE, MSE, coverage50, coverage95 in percent, WIS, WIS_log).
%                   With 'Aggregate' 1 the MAE, MSE, coverage95 and WIS columns are
%                   the values of the 'performance-forecasting' files written by
%                   the runners (computeforecastperformance and computeWIS).
%     S.summary     struct, means over all scored horizons: nScored, MAE, MSE,
%                   RMSE, coverage50, coverage95, WIS, WIS_log and, with a
%                   baseline, WIS_baseline and relWIS = mean(WIS)/mean(WIS_baseline)
%
%   Definitions (identical to computeWIS and computeforecastperformance)
%     point forecast  the median of the draws
%     interval        quantiles of the draws, set to zero when negative
%     WIS             (1/(K + 1/2)) * ( 0.5*|y - median| + sum_k (alpha_k/2) * IS_k ),
%                     alpha = [0.02 0.05 0.1:0.1:0.9], K = 11
%   A horizon without an observation gets NaN scores and is left out of the
%   summary; the forecast-hub convention of scoring each horizon separately and
%   averaging the scores is S.perHorizon followed by S.summary.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

p = inputParser;
p.addParameter('Aggregate', 1);
p.addParameter('BaselineDraws', []);
p.parse(varargin{:});
k = p.Results.Aggregate;
baseDraws = p.Results.BaselineDraws;

if ~isnumeric(draws) || ~ismatrix(draws) || isempty(draws) || ~isreal(draws)
    error('quantdiffForecastScores:draws', 'draws must be a real numeric matrix, one column per trajectory.');
end
if ~isscalar(nCalib) || nCalib < 0 || nCalib ~= fix(nCalib) || ~isscalar(horizon) || horizon < 1 || horizon ~= fix(horizon)
    error('quantdiffForecastScores:sizes', 'nCalib must be a nonnegative integer and horizon a positive integer.');
end
if size(draws, 1) < nCalib + horizon
    error('quantdiffForecastScores:rows', ...
        'draws has %d rows; %d calibration rows + %d forecast steps are needed.', size(draws, 1), nCalib, horizon);
end
if ~isscalar(k) || k < 1 || k ~= fix(k)
    error('quantdiffForecastScores:aggregate', '''Aggregate'' must be a positive integer.');
end
if size(datalatest, 2) < 2
    error('quantdiffForecastScores:data', 'datalatest must be [time value].');
end
nB = floor(horizon / k);                       % horizons (blocks) that are scored
if nB < 1
    error('quantdiffForecastScores:aggregate', ...
        '''Aggregate'' = %d is longer than the horizon (%d steps).', k, horizon);
end
if ~isempty(baseDraws) && size(baseDraws, 1) < nCalib + horizon
    error('quantdiffForecastScores:baseline', ...
        'BaselineDraws has %d rows; %d are needed.', size(baseDraws, 1), nCalib + horizon);
end

% Forecast-period draws, summed over blocks of k steps.
D = localAggregate(draws(nCalib + (1:nB * k), :), k, nB);

% Observations summed over the same blocks (NaN where the data end).
nObs = size(datalatest, 1);
y = nan(nB, 1); tFrom = nan(nB, 1); tTo = nan(nB, 1);
stepFrom = ((0:nB - 1) * k + 1)';
stepTo = stepFrom + k - 1;
for b = 1:nB
    rows = nCalib + (stepFrom(b):stepTo(b));
    if rows(1) <= nObs, tFrom(b) = datalatest(rows(1), 1); end
    if rows(end) <= nObs
        tTo(b) = datalatest(rows(end), 1);
        y(b) = sum(datalatest(rows, 2));
    end
end

[wis, med, c50, c95] = localScore(D, y);
wisLog = localScore(log1p(max(D, 0)), log1p(max(y, 0)));

absErr = abs(y - med);
sqErr = (y - med) .^ 2;
scored = isfinite(y);
c50(~scored) = NaN; c95(~scored) = NaN;
wisLog(~scored) = NaN;                         % max(NaN, 0) is 0 in MATLAB, so set these explicitly

S.perHorizon = table((1:nB)', stepFrom, stepTo, tFrom, tTo, y, med, absErr, sqErr, c50, c95, wis, wisLog, ...
    'VariableNames', {'horizon', 'step_from', 'step_to', 'time_from', 'time_to', 'observed', 'median', ...
    'abs_error', 'sq_error', 'cover50', 'cover95', 'WIS', 'WIS_log'});

hasBase = ~isempty(baseDraws);
if hasBase
    wisBase = localScore(localAggregate(baseDraws(nCalib + (1:nB * k), :), k, nB), y);
    S.perHorizon.WIS_baseline = wisBase;
    S.perHorizon.relWIS = wis ./ wisBase;
end

% Running means over horizons 1..h (the convention of the runners' files).
h = (1:nB)';
S.cumulative = table(h, cumsum(absErr) ./ h, cumsum(sqErr) ./ h, 100 * cumsum(c50) ./ h, 100 * cumsum(c95) ./ h, ...
    cumsum(wis) ./ h, cumsum(wisLog) ./ h, ...
    'VariableNames', {'horizon', 'MAE', 'MSE', 'coverage50', 'coverage95', 'WIS', 'WIS_log'});

% Means over the scored horizons.
S.summary = struct('nScored', nnz(scored), 'MAE', mean(absErr(scored)), 'MSE', mean(sqErr(scored)), ...
    'RMSE', sqrt(mean(sqErr(scored))), 'coverage50', 100 * mean(c50(scored)), 'coverage95', 100 * mean(c95(scored)), ...
    'WIS', mean(wis(scored)), 'WIS_log', mean(wisLog(scored)));
if hasBase
    both = scored & isfinite(wisBase);
    S.summary.WIS_baseline = mean(wisBase(both));
    S.summary.relWIS = mean(wis(both)) / mean(wisBase(both));
end
end


% <============================================================================>
% < Sum every trajectory over blocks of k consecutive rows                    ==>
% <============================================================================>
function A = localAggregate(F, k, nB)
M = size(F, 2);
if k == 1
    A = F;
else
    A = reshape(sum(reshape(F, k, nB, M), 1), nB, M);
end
end


% <============================================================================>
% < WIS, median and interval coverage of each row of D against y              ==>
% <============================================================================>
function [wis, med, c50, c95] = localScore(D, y)
alphas = [0.02 0.05 0.1:0.1:0.9];
K = numel(alphas);
n = size(D, 1);
wis = nan(n, 1); med = nan(n, 1); c50 = nan(n, 1); c95 = nan(n, 1);
levels = [alphas / 2, 1 - alphas / 2, 0.5, 0.25, 0.75, 0.025, 0.975];
for t = 1:n
    q = quantile(D(t, :), levels);
    L = max(q(1:K), 0);
    U = max(q(K + 1:2 * K), 0);
    med(t) = q(2 * K + 1);
    IS = (U - L) + (2 ./ alphas) .* (L - y(t)) .* (y(t) < L) + (2 ./ alphas) .* (y(t) - U) .* (y(t) > U);
    wis(t) = (0.5 * abs(y(t) - med(t)) + sum((alphas / 2) .* IS)) / (K + 0.5);
    c50(t) = double(y(t) >= max(q(2 * K + 2), 0) && y(t) <= max(q(2 * K + 3), 0));
    c95(t) = double(y(t) >= max(q(2 * K + 4), 0) && y(t) <= max(q(2 * K + 5), 0));
end
end
