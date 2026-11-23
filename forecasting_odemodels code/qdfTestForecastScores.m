function results = qdfTestForecastScores()
%QDFTESTFORECASTSCORES  Checks of the per-horizon forecast scores (quantdiffForecastScores).
%
%   qdfTestForecastScores            run the nine checks and print a verdict each
%   results = qdfTestForecastScores  also return {name, verdict, detail}
%
%   Run it from the toolbox code folder. No model is fitted, so it takes a few
%   seconds. It needs the Statistics and Machine Learning Toolbox (quantile).
%   The checks use a fixed matrix of predictive draws (5 calibration rows,
%   6 forecast steps, 201 draws) and compare with an independent Python/SciPy
%   calculation on the same matrix.
%
%   Checks
%     1  the running means of the per-horizon scores equal the existing scores
%        (computeforecastperformance and computeWIS): MAE, MSE, 95% coverage, WIS
%     2  per-horizon WIS and log-scale WIS agree with Python
%     3  50% and 95% coverage per horizon, and the summary
%     4  totals over blocks of 3 steps ('Aggregate', 3) agree with Python
%     5  the total over the whole horizon equals scoring the summed trajectories,
%        and is not the sum of the daily errors
%     6  relative WIS: 1 against itself; agrees with Python against a second
%        forecast; the summary is the ratio of the mean scores
%     7  a perfect forecast scores 0 with full coverage
%     8  horizons beyond the data get NaN and are left out of the summary
%     9  the flat baseline can be scored, and quantdiffCompareForecasts stops with
%        a clear message for a forecast that has not been run
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

global method1 %#ok<GVMIS>
savedMethod = method1;
restore = onCleanup(@() localRestore(savedMethod)); %#ok<NASGU>

z = linspace(-2.5, 2.5, 201);
mu = [10 12 15 20 26 36 45 3 70 90 115]';
draws = mu + 1.5 * sqrt(mu) * z;                       % 11 x 201
obs = [10 11 16 19 27 40 44 0 60 130 100]';
datalatest = [(0:10)' obs];
nCalib = 5; H = 6;
data1 = datalatest(1:nCalib, :);

refWIS  = [3.7014498913 3.81834247523 1.49707346475 6.10643287028 25.3475931643 8.49112849215]';
refLog  = [0.103663461034 0.0859728985552 0.595357849149 0.0895407128021 0.23058620871 0.0759060403564]';
refAgg3 = [8.1258711591 17.1010805564]';
refAgg6 = 24.8591296959;
refBase = [3.6442033913 3.23532411186 5.25987900115 3.80994454666 41.2780745314 4.95806486139]';
refRel  = 0.787354413522;
tol = 1e-8;

results = cell(0, 3);
S = quantdiffForecastScores(draws, datalatest, nCalib, H);
P = S.perHorizon;

% 1 ---------------------------------------------------------------------------
[~, ~, ~, ~, ~, ~, ~, MSEF, MAEF, PIF] = evalc('computeforecastperformance(data1, datalatest, draws, draws, H)');
[~, ~, WISF] = evalc('computeWIS(data1, datalatest, draws, H)');
d1 = [max(abs(S.cumulative.MAE - MAEF(:, 2))), max(abs(S.cumulative.MSE - MSEF(:, 2))), ...
      max(abs(S.cumulative.coverage95 - PIF(:, 2))), max(abs(S.cumulative.WIS - WISF(:, 2)))];
results = localAdd(results, '1 running means = existing scores', all(d1 < 1e-9), sprintf('largest differences: %s', mat2str(d1, 3)));

% 2 ---------------------------------------------------------------------------
d2 = [max(abs(P.WIS - refWIS) ./ refWIS), max(abs(P.WIS_log - refLog) ./ refLog)];
results = localAdd(results, '2 per-horizon WIS vs Python', all(d2 < tol), ...
    sprintf('WIS = %s; largest relative differences %s', mat2str(P.WIS', 6), mat2str(d2, 3)));

% 3 ---------------------------------------------------------------------------
ok3 = isequal(P.cover50', [1 1 1 1 0 1]) && isequal(P.cover95', [1 1 1 1 0 1]) && isequal(P.observed', [40 44 0 60 130 100]) && ...
    max(abs(P.median' - [36 45 3 70 90 115])) < 1e-9 && S.summary.nScored == 6 && ...
    abs(S.summary.WIS - 8.16033672633) < 1e-8 && abs(S.summary.coverage95 - 500 / 6) < 1e-9 && abs(S.summary.MAE - 73 / 6) < 1e-9;
results = localAdd(results, '3 coverage and summary', ok3, ...
    sprintf('cover95 = %s; mean WIS = %.8f; MAE = %.4f', mat2str(P.cover95'), S.summary.WIS, S.summary.MAE));

% 4 ---------------------------------------------------------------------------
S3 = quantdiffForecastScores(draws, datalatest, nCalib, H, 'Aggregate', 3);
ok4 = height(S3.perHorizon) == 2 && isequal(S3.perHorizon.observed', [84 290]) && max(abs(S3.perHorizon.median' - [84 275])) < 1e-9 && ...
    isequal(S3.perHorizon.step_from', [1 4]) && isequal(S3.perHorizon.step_to', [3 6]) && ...
    max(abs(S3.perHorizon.WIS - refAgg3) ./ refAgg3) < tol;
results = localAdd(results, '4 totals over blocks of 3', ok4, ...
    sprintf('observed = %s; WIS = %s', mat2str(S3.perHorizon.observed'), mat2str(S3.perHorizon.WIS', 8)));

% 5 ---------------------------------------------------------------------------
S6 = quantdiffForecastScores(draws, datalatest, nCalib, H, 'Aggregate', 6);
Sd = quantdiffForecastScores(sum(draws(nCalib + 1:end, :), 1), [0 374], 0, 1);
S4 = quantdiffForecastScores(draws, datalatest, nCalib, H, 'Aggregate', 4);     % incomplete second block is dropped
threw = false;
try
    quantdiffForecastScores(draws, datalatest, nCalib, H, 'Aggregate', 7);
catch ME
    threw = strcmp(ME.identifier, 'quantdiffForecastScores:aggregate');
end
ok5 = S6.perHorizon.observed == 374 && abs(S6.perHorizon.WIS - refAgg6) < tol * refAgg6 && ...
    abs(Sd.perHorizon.WIS - S6.perHorizon.WIS) < 1e-9 && abs(S6.perHorizon.abs_error - 15) < 1e-9 && abs(sum(P.abs_error) - 73) < 1e-9 && ...
    height(S4.perHorizon) == 1 && S4.perHorizon.observed == 144 && threw;
results = localAdd(results, '5 total over the horizon', ok5, ...
    sprintf('WIS of the total = %.8f; error of the total = %g, sum of daily errors = %g', S6.perHorizon.WIS, S6.perHorizon.abs_error, sum(P.abs_error)));

% 6 ---------------------------------------------------------------------------
Ssame = quantdiffForecastScores(draws, datalatest, nCalib, H, 'BaselineDraws', draws);
Sb = quantdiffForecastScores(draws, datalatest, nCalib, H, 'BaselineDraws', 0.8 * draws + 5);
ok6 = all(Ssame.perHorizon.relWIS == 1) && Ssame.summary.relWIS == 1 && ...
    max(abs(Sb.perHorizon.WIS_baseline - refBase) ./ refBase) < tol && abs(Sb.summary.relWIS - refRel) < tol && ...
    abs(Sb.summary.relWIS - mean(P.WIS) / mean(Sb.perHorizon.WIS_baseline)) < 1e-12;
results = localAdd(results, '6 relative WIS', ok6, sprintf('relative WIS = %.9f (Python %.9f)', Sb.summary.relWIS, refRel));

% 7 ---------------------------------------------------------------------------
yP = [5 0 120]';
Sp = quantdiffForecastScores(repmat(yP, 1, 50), [(0:2)' yP], 0, 3);
ok7 = all(Sp.perHorizon.WIS == 0) && all(Sp.perHorizon.WIS_log == 0) && Sp.summary.coverage50 == 100 && ...
    Sp.summary.coverage95 == 100 && Sp.summary.MAE == 0;
results = localAdd(results, '7 perfect forecast', ok7, sprintf('WIS = %s, coverage95 = %g', mat2str(Sp.perHorizon.WIS'), Sp.summary.coverage95));

% 8 ---------------------------------------------------------------------------
Sshort = quantdiffForecastScores(draws, datalatest(1:nCalib + 4, :), nCalib, H);
T8 = Sshort.perHorizon;
ok8 = all(isfinite(T8.WIS(1:4))) && all(isnan(T8.WIS(5:6))) && all(isnan(T8.WIS_log(5:6))) && all(isnan(T8.cover95(5:6))) && ...
    Sshort.summary.nScored == 4 && abs(Sshort.summary.WIS - mean(refWIS(1:4))) < 1e-8;
results = localAdd(results, '8 horizons beyond the data', ok8, sprintf('scored %d of %d; mean WIS = %.8f', Sshort.summary.nScored, H, Sshort.summary.WIS));

% 9 ---------------------------------------------------------------------------
baseDraws = baselineForecast(data1, H, 'flat', 400, 1);
Sf = quantdiffForecastScores(draws, datalatest, nCalib, H, 'BaselineDraws', baseDraws);
medBase = median(baseDraws(nCalib + 1:end, :), 2);
missingCaught = false;
try
    evalc('quantdiffCompareForecasts(@options_forecast_SEIR_flu1918_dist1_3, [46 7], 3, ''Save'', false)');
catch ME
    missingCaught = strcmp(ME.identifier, 'quantdiffCompareForecasts:missing');
end
badCaught = false;
try
    quantdiffCompareForecasts(@options_forecast_SEIR_flu1918_dist1_3, [1 17 3], 10);
catch ME
    badCaught = strcmp(ME.identifier, 'quantdiffCompareForecasts:windows');
end
ok9 = all(isfinite(Sf.perHorizon.WIS_baseline)) && all(Sf.perHorizon.WIS_baseline > 0) && ...
    max(abs(medBase - obs(nCalib))) < 1e-9 && isfinite(Sf.summary.relWIS) && missingCaught && badCaught;
results = localAdd(results, '9 baseline and comparison checks', ok9, ...
    sprintf('baseline median = %g at every step (last calibration value %g); relative WIS = %.4f; errors caught: %d %d', ...
    medBase(1), obs(nCalib), Sf.summary.relWIS, missingCaught, badCaught));

fprintf('\n%d of %d checks passed.\n', sum(strcmp(results(:, 2), 'PASS')), size(results, 1));
end


function results = localAdd(results, name, ok, detail)
if ok, verdict = 'PASS'; else, verdict = 'FAIL'; end
fprintf('%-36s %s   %s\n', name, verdict, detail);
results(end+1, :) = {name, verdict, detail};
end

function localRestore(savedMethod)
global method1 %#ok<GVMIS>
method1 = savedMethod;
end
