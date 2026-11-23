function R = quantdiffCompareForecasts(optionsHandles, windows, horizon, varargin)
%QUANTDIFFCOMPAREFORECASTS  Score several models on the same targets, by horizon and over several forecast origins.
%
%   R = quantdiffCompareForecasts({@options_forecast_SEIR_flu1918_dist1_3, @options_forecast_SEIR_flu1918_IC}, ...
%                                 [1 14; 1 17; 1 20], 10)
%   R = quantdiffCompareForecasts(..., 'Labels', {'quick start', 'I_0 estimated'}, 'Run', true)
%   R = quantdiffCompareForecasts(..., 'Aggregate', 7, 'Baseline', 'trend')
%
%   optionsHandles  cell of forecast-options function handles (the ones passed
%                   to Run_Forecasting_ODEModel), or a single handle
%   windows         one row per forecast origin: [tstart windowsize]. The forecast
%                   is made from calibration rows tstart .. tstart+windowsize-1
%                   and scored on the rows that follow. [1 14; 1 17; 1 20] is an
%                   expanding window; [1 17; 2 17; 3 17] a rolling one.
%   horizon         forecast horizon the models were run with
%
%   Options (name/value)
%     'Labels'     cell of names for the models (default: the options function names)
%     'Tend'       the '-tend-' value in the result file names: tend1 of the run that
%                  produced them. Scalar or one value per window. Default: equal to
%                  tstart, i.e. Run_Forecasting_ODEModel(opts, tstart, tstart, W, horizon).
%     'Run'        true -> a forecast whose result file is missing is run first, with
%                  rng('Seed') set before each run. Default false: the function stops
%                  and prints the command to run.
%     'Seed'       seed for those runs and for the baseline draws (default 1)
%     'Aggregate'  k: score totals over blocks of k forecast steps (default 1)
%     'Baseline'   'flat' (default), 'trend' or 'none' (see baselineForecast)
%     'FitIndex'   which fitted series to score (default 1)
%     'Save'       true (default) -> write three CSV files to <toolbox>/output
%     'Name'       tag used in those file names (default 'comparison')
%
%   The result files are found from the options functions with
%   qdf_forecast_matname, so the model, method, window and horizon of every file
%   are exactly the ones asked for; nothing is matched by pattern. All models must
%   have been fitted to the same observed series.
%
%   Output R (struct)
%     R.perHorizon  table, one row per model, window and horizon, each scored on
%                   its own (quantdiffForecastScores): model, tstart, windowsize,
%                   origin_time, horizon, observed, median, abs_error, sq_error,
%                   cover50, cover95, WIS, WIS_log, WIS_baseline, relWIS
%     R.byHorizon   table, one row per model and horizon: means over the windows
%                   (MAE, RMSE, coverage50 and coverage95 in percent, WIS, WIS_log,
%                   WIS_baseline, relWIS = mean WIS / mean baseline WIS, nWindows)
%     R.summary     table, one row per model: the same means over all windows and
%                   horizons (nTargets = number of targets averaged)
%     R.files       result files read (windows x models)
%   The means use only the targets that every model has a score for, so the
%   models are always compared on the same observations. The mean of the WIS
%   over horizons and origins, divided by the baseline's, is the relative WIS of
%   the forecast hubs.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

global method1 %#ok<GVMIS>  set as a side effect inside the options functions

if isa(optionsHandles, 'function_handle'), optionsHandles = {optionsHandles}; end
if ~iscell(optionsHandles) || isempty(optionsHandles) || ~all(cellfun(@(h) isa(h, 'function_handle'), optionsHandles))
    error('quantdiffCompareForecasts:options', 'optionsHandles must be a function handle or a nonempty cell array of function handles.');
end
K = numel(optionsHandles);
if ~isnumeric(windows) || isempty(windows) || size(windows, 2) ~= 2 || any(windows(:) < 1) || any(windows(:) ~= fix(windows(:)))
    error('quantdiffCompareForecasts:windows', 'windows must have one row per forecast origin: [tstart windowsize], positive integers.');
end
if ~isscalar(horizon) || horizon < 1 || horizon ~= fix(horizon)
    error('quantdiffCompareForecasts:horizon', 'horizon must be a positive integer.');
end
nW = size(windows, 1);

p = inputParser;
p.addParameter('Labels', {});
p.addParameter('Tend', []);
p.addParameter('Run', false);
p.addParameter('Seed', 1);
p.addParameter('Aggregate', 1);
p.addParameter('Baseline', 'flat');
p.addParameter('FitIndex', 1);
p.addParameter('Save', true);
p.addParameter('Name', 'comparison');
p.parse(varargin{:});
labels = p.Results.Labels; tendIn = p.Results.Tend; doRun = p.Results.Run; seed = p.Results.Seed;
agg = p.Results.Aggregate; baseType = lower(char(p.Results.Baseline)); fitIdx = p.Results.FitIndex;
doSave = p.Results.Save; tag = char(p.Results.Name);

if isempty(labels)
    labels = cellfun(@func2str, optionsHandles, 'UniformOutput', false);
end
labels = cellstr(labels);
if numel(labels) ~= K
    error('quantdiffCompareForecasts:labels', '''Labels'' must have one entry per options handle (%d).', K);
end
if isempty(tendIn)
    tend = windows(:, 1);
elseif isscalar(tendIn)
    tend = repmat(tendIn, nW, 1);
else
    tend = tendIn(:);
    if numel(tend) ~= nW
        error('quantdiffCompareForecasts:tend', '''Tend'' must be a scalar or have one value per window (%d).', nW);
    end
end
if ~any(strcmp(baseType, {'flat', 'trend', 'none'}))
    error('quantdiffCompareForecasts:baseline', '''Baseline'' must be ''flat'', ''trend'' or ''none''.');
end

basePath = fileparts(mfilename('fullpath'));
outputDir = fullfile(basePath, 'output');

nH = floor(horizon / agg);                       % horizons scored per window
if nH < 1
    error('quantdiffCompareForecasts:aggregate', '''Aggregate'' = %d is longer than the horizon (%d steps).', agg, horizon);
end

% Per-target scores: windows x horizons x models.
AE = nan(nW, nH, K); SE = AE; C50 = AE; C95 = AE; WIS = AE; WISL = AE; MED = AE;
OBS = nan(nW, nH); WISB = nan(nW, nH); ORIGIN = nan(nW, 1);
R.files = cell(nW, K);

for w = 1:nW
    tstart = windows(w, 1); W = windows(w, 2);
    obsRef = [];
    for k = 1:K
        [cad, ~, ~, dist1k, ~, ~, modelk, paramsk, varsk] = optionsHandles{k}();
        if ~isempty(method1) && method1 > 0, dist1k = method1; end
        fname = qdf_forecast_matname(cad, modelk.name, varsk.fit_index(fitIdx), paramsk.fixI0, ...
            method1, dist1k, tstart, tend(w), W, horizon);
        fpath = fullfile(outputDir, fname);
        if ~isfile(fpath) && doRun && tend(w) == tstart
            fprintf('quantdiffCompareForecasts: running %s for rows %d-%d, horizon %d ...\n', labels{k}, tstart, tstart + W - 1, horizon);
            rng(seed, 'twister');
            Run_Forecasting_ODEModel(optionsHandles{k}, tstart, tstart, W, horizon);
        end
        if ~isfile(fpath)
            error('quantdiffCompareForecasts:missing', ...
                ['No result file for %s, calibration rows %d-%d, horizon %d:\n  %s\n', ...
                 'Run  Run_Forecasting_ODEModel(@%s, %d, %d, %d, %d)  first, or call this function with ''Run'', true.'], ...
                labels{k}, tstart, tstart + W - 1, horizon, fpath, func2str(optionsHandles{k}), tstart, tend(w), W, horizon);
        end
        F = load(fpath, 'forecast2', 'datalatest', 'data1');
        if ~all(isfield(F, {'forecast2', 'datalatest', 'data1'}))
            error('quantdiffCompareForecasts:vars', '%s lacks forecast2, datalatest or data1 (written by an older runner?).', fname);
        end
        R.files{w, k} = fname;
        nCalib = size(F.data1, 1);
        obs = F.datalatest(:, [1 fitIdx + 1]);
        calib = F.data1(:, [1 fitIdx + 1]);

        % All models must forecast the same observed series from the same origin.
        nCmp = min(size(obs, 1), nCalib + horizon);
        if k == 1
            obsRef = obs(1:nCmp, :);
            ORIGIN(w) = obs(nCalib, 1);
            baseDraws = [];
            if ~strcmp(baseType, 'none')
                baseDraws = baselineForecast(calib, horizon, baseType, size(F.forecast2, 2), seed);
            end
        elseif nCalib ~= W || ~isequal(obs(1:nCmp, :), obsRef)
            error('quantdiffCompareForecasts:data', ...
                '%s was fitted to different data than %s for calibration rows %d-%d; the models must share the observed series.', ...
                labels{k}, labels{1}, tstart, tstart + W - 1);
        end

        if isempty(baseDraws)
            S = quantdiffForecastScores(F.forecast2, obs, nCalib, horizon, 'Aggregate', agg);
        else
            S = quantdiffForecastScores(F.forecast2, obs, nCalib, horizon, 'Aggregate', agg, 'BaselineDraws', baseDraws);
            WISB(w, :) = S.perHorizon.WIS_baseline';
        end
        T = S.perHorizon;
        OBS(w, :) = T.observed';
        MED(w, :, k) = T.median';   AE(w, :, k) = T.abs_error';  SE(w, :, k) = T.sq_error';
        C50(w, :, k) = T.cover50';  C95(w, :, k) = T.cover95';
        WIS(w, :, k) = T.WIS';      WISL(w, :, k) = T.WIS_log';
    end
end

% ---- long table: one row per model, window and horizon
nRows = K * nW * nH;
cModel = cell(nRows, 1); M = nan(nRows, 13); r = 0;
for k = 1:K
    for w = 1:nW
        for h = 1:nH
            r = r + 1;
            cModel{r} = labels{k};
            M(r, :) = [windows(w, 1) windows(w, 2) ORIGIN(w) h OBS(w, h) MED(w, h, k) AE(w, h, k) SE(w, h, k) ...
                C50(w, h, k) C95(w, h, k) WIS(w, h, k) WISL(w, h, k) WISB(w, h)];
        end
    end
end
R.perHorizon = [table(cModel, 'VariableNames', {'model'}), array2table(M, 'VariableNames', ...
    {'tstart', 'windowsize', 'origin_time', 'horizon', 'observed', 'median', 'abs_error', 'sq_error', ...
     'cover50', 'cover95', 'WIS', 'WIS_log', 'WIS_baseline'})];
R.perHorizon.relWIS = R.perHorizon.WIS ./ R.perHorizon.WIS_baseline;

% ---- means over the targets every model has a score for
common = all(isfinite(WIS), 3);                       % windows x horizons
useBase = ~strcmp(baseType, 'none');
if useBase, common = common & isfinite(WISB); end

cModel = cell(K * nH, 1); B = nan(K * nH, 10); r = 0;
for k = 1:K
    for h = 1:nH
        r = r + 1;
        cModel{r} = labels{k};
        B(r, :) = [h, localMeans(AE(:, h, k), SE(:, h, k), C50(:, h, k), C95(:, h, k), WIS(:, h, k), WISL(:, h, k), WISB(:, h), common(:, h))];
    end
end
namesB = {'horizon', 'MAE', 'RMSE', 'coverage50', 'coverage95', 'WIS', 'WIS_log', 'WIS_baseline', 'relWIS', 'nWindows'};
R.byHorizon = [table(cModel, 'VariableNames', {'model'}), array2table(B, 'VariableNames', namesB)];

Sm = nan(K, 9);
for k = 1:K
    Sm(k, :) = localMeans(reshape(AE(:, :, k), [], 1), reshape(SE(:, :, k), [], 1), reshape(C50(:, :, k), [], 1), ...
        reshape(C95(:, :, k), [], 1), reshape(WIS(:, :, k), [], 1), reshape(WISL(:, :, k), [], 1), WISB(:), common(:));
end
R.summary = [table(labels(:), 'VariableNames', {'model'}), array2table(Sm, 'VariableNames', ...
    {'MAE', 'RMSE', 'coverage50', 'coverage95', 'WIS', 'WIS_log', 'WIS_baseline', 'relWIS', 'nTargets'})];
R.windows = windows; R.horizon = horizon; R.aggregate = agg; R.baseline = baseType;

if doSave
    if ~isfolder(outputDir), mkdir(outputDir); end
    stem = fullfile(outputDir, ['scores-comparison-' regexprep(tag, '[\s/\\:*?"<>|]+', '_')]);
    writetable(R.perHorizon, [stem '-perhorizon.csv']);
    writetable(R.byHorizon, [stem '-byhorizon.csv']);
    writetable(R.summary, [stem '-summary.csv']);
    R.savedTo = stem;
end

fprintf('\nForecast comparison: %d window(s), horizon %d', nW, horizon);
if agg > 1, fprintf(', totals over blocks of %d steps', agg); end
fprintf(', baseline: %s. Means over the %d target(s) scored for every model.\n', baseType, nnz(common));
disp(R.summary)
end


% <============================================================================>
% < Means over the selected targets                                           ==>
% <============================================================================>
function v = localMeans(ae, se, c50, c95, wis, wisl, wisb, use)
% [MAE RMSE coverage50 coverage95 WIS WIS_log WIS_baseline relWIS n]
use = logical(use(:));
n = nnz(use);
if n == 0
    v = [nan(1, 8) 0];
    return
end
wb = mean(wisb(use));
v = [mean(ae(use)), sqrt(mean(se(use))), 100 * mean(c50(use)), 100 * mean(c95(use)), ...
     mean(wis(use)), mean(wisl(use)), wb, mean(wis(use)) / wb, n];
end
