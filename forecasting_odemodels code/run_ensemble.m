function R = run_ensemble(optionsHandles, modelNames, weight_type1, calibrationperiods1, forecastingperiod1, name_ensemble, varargin)
%RUN_ENSEMBLE  Weighted linear-pool ensemble of K previously run QuantDiffForecast models.
%
%   R = run_ensemble({@options_forecast_GLM_x, @options_forecast_RICH_x}, {'GLM','RICH'}, ...
%                    weight_type1, calibrationperiods1, forecastingperiod1, 'GLM-RICH-Ensemble')
%   R = run_ensemble(..., 'tstart', 1, 'tend', 1, 'Nens', [], 'Seed', 1, 'Plot', true)
%
%   Replaces the six temp_run_ENS*.m scripts. For every calibration window it
%   (1) loads each model's 'Forecast-ODEModel-...' result written by
%       Run_Forecasting_ODEModel (file name from qdf_forecast_matname, so the
%       reader cannot drift from the writer), reading only the variables it
%       needs instead of dumping the saved workspace into this function;
%   (2) computes the weights with ensembleWeights (equal, Akaike, inverse-WIS or
%       softmax-WIS; the better model always gets the larger weight);
%   (3) pools the models' predictive draws with the exact empirical linear pool
%       getensemble_linearpool_from_curves and scores the ensemble from the pool's
%       EXACT quantiles (computeScoresFromQuantiles), so its scores carry no
%       sampling noise and do not depend on 'Seed'. Members and baselines are
%       scored through the same function from the exact quantiles of their own
%       draws. The pool's trajectory sample curvesEns only draws the figure;
%   (4) writes the ensemble Forecast / performance CSVs to <toolbox>/output and
%       returns everything in R.
%
%   Inputs
%     optionsHandles      cell of K options function handles (the same ones used
%                         to run Run_Forecasting_ODEModel for each model)
%     modelNames          cell of K display names
%     weight_type1        -1 equal | 0/1 Akaike | 2 inverse WIS | 3 softmax WIS
%     calibrationperiods1 vector of window sizes already run
%     forecastingperiod1  forecast horizon already run
%     name_ensemble       label used in output file names
%   Options (name/value)
%     'tstart','tend'     the '-tstart-' token of the result file is the index of ONE
%                         rolling window (the runner writes one file per window), so
%                         'tstart' selects that window and 'tend' must equal the tend1 the
%                         models were run with (defaults 1, 1)
%     'Nens'              number of pooled trajectories drawn for the figure, and of
%                         baseline draws (default: M of the members)
%     'Seed'              seed of the trajectory sample and of the baselines (default 1);
%                         the ensemble's and the members' scores do not depend on it
%     'Plot'              true/false (default true)
%     'FitIndex'          which fitted variable's file to read (default 1)
%
%   Output R (struct): calibrationperiods, weights (nW x K), rule, members
%   (names), performanceC_models / performanceF_models (nW x 4 x K: MAE MSE
%   Coverage WIS), performanceC_ens / performanceF_ens (nW x 4; the latter at the
%   full horizon), performanceF_ens_byhorizon (1 x nW cell; rows [h MAE MSE Coverage
%   WIS WIS_baseline_flat WIS_baseline_trend relWIS_flat relWIS_trend] for
%   h = 1..horizon, each averaged over horizons 1..h like the members' files),
%   WIS_baselines (nW x 2: hub flat, log-linear trend; see baselineForecast),
%   relWIS_ens (nW x 2), relWIS_models (nW x K x 2) at the full horizon,
%   AICcs (nW x K), files (the .mat files read).
%
%   Writer-side note: the runners build the same file name inline
%   (the 'Forecast-ODEModel-' strcat near the per-window save in
%   Run_Forecasting_ODEModel.m and Run_Fit_ODEModel.m). Replacing that
%   strcat with qdf_forecast_matname(cadfilename1, model.name, vars.fit_index(j),
%   params.fixI0, method1, dist1, i, tend1, windowsize1, forecastingperiod)
%   keeps writer and reader in one place; until then the helper reproduces the
%   current convention exactly.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

global method1 %#ok<GVMIS>  set as a side effect inside the options functions

p = inputParser;
p.addParameter('tstart', 1);
p.addParameter('tend', 1);
p.addParameter('Nens', []);
p.addParameter('Seed', 1);
p.addParameter('Plot', true);
p.addParameter('FitIndex', 1);
p.parse(varargin{:});
tstart1 = p.Results.tstart; tend1 = p.Results.tend;
Nens = p.Results.Nens; seed = p.Results.Seed; doPlot = p.Results.Plot; fitIdx = p.Results.FitIndex;

if ~iscell(optionsHandles) || isempty(optionsHandles) || ~all(cellfun(@(h) isa(h, 'function_handle'), optionsHandles))
    error('run_ensemble:options', 'optionsHandles must be a nonempty cell array of function handles.');
end
K = numel(optionsHandles);
if numel(modelNames) ~= K
    error('run_ensemble:names', 'modelNames must have one entry per options handle (%d).', K);
end

basePath = fileparts(mfilename('fullpath'));
outputDir = fullfile(basePath, 'output');
if ~isfolder(outputDir), mkdir(outputDir); end

nW = numel(calibrationperiods1);
R = struct('calibrationperiods', calibrationperiods1(:)', 'members', {modelNames}, 'weight_type1', weight_type1, ...
    'rule', '', 'weights', nan(nW, K), 'AICcs', nan(nW, K), ...
    'performanceC_models', nan(nW, 4, K), 'performanceF_models', nan(nW, 4, K), ...
    'performanceC_ens', nan(nW, 4), 'performanceF_ens', nan(nW, 4), 'performanceF_ens_byhorizon', {cell(1, nW)}, ...
    'WIS_baselines', nan(nW, 2), 'relWIS_ens', nan(nW, 2), 'relWIS_models', nan(nW, K, 2), 'files', {cell(nW, K)});

needVars = {'forecast2', 'AICcs', 'MAECSS', 'MSECSS', 'PICSS', 'WISCSS', 'MAEFSS', 'MSEFSS', 'PIFSS', 'WISFSS', ...
            'timevect1', 'timevect2', 'timevect_all', 'data_all', 'data1', 'datalatest', 'forecastingperiod'};

for w = 1:nW
    windowsize1 = calibrationperiods1(w);
    curves = cell(1, K); AICc = nan(1, K); WISC = nan(1, K); S = cell(1, K);
    caddisease1 = ''; datatype1 = ''; dist1_1 = NaN;

    for k = 1:K
        [cadfilename1_INP, caddisease_INP, datatype_INP, dist1_INP, ~, ~, model_INP, params_INP, vars_INP] = optionsHandles{k}();
        dist1_k = dist1_INP; if ~isempty(method1) && method1 > 0, dist1_k = method1; end
        fname = qdf_forecast_matname(cadfilename1_INP, model_INP.name, vars_INP.fit_index(fitIdx), params_INP.fixI0, ...
                                     method1, dist1_k, tstart1, tend1, windowsize1, forecastingperiod1);
        fpath = fullfile(outputDir, fname);
        if ~isfile(fpath)
            error('run_ensemble:missing', ['Result file for model %d (%s), window %d not found:\n  %s\n', ...
                'Run Run_Forecasting_ODEModel(%s, %d, %d, %d, %d) first.'], k, modelNames{k}, windowsize1, fpath, ...
                func2str(optionsHandles{k}), tstart1, tend1, windowsize1, forecastingperiod1);
        end
        S{k} = load(fpath, needVars{:});                  % named variables only: no workspace dump
        missing = needVars(~isfield(S{k}, needVars));
        if ~isempty(missing)
            error('run_ensemble:vars', '%s lacks variables %s (written by an older runner?).', fname, strjoin(missing, ', '));
        end
        curves{k} = max(S{k}.forecast2, 0);                % the per-variable slice saved with this file (forecast_model12 stacks ALL fitted variables)
        AICc(k) = S{k}.AICcs(end, 2);                      % AICcs accumulates over rolling windows; the last row is this window's
        WISC(k) = S{k}.WISCSS(end);
        R.performanceC_models(w, :, k) = [S{k}.MAECSS(end) S{k}.MSECSS(end) S{k}.PICSS(end) S{k}.WISCSS(end)];
        if ~isempty(S{k}.MAEFSS)
            R.performanceF_models(w, :, k) = [S{k}.MAEFSS(end) S{k}.MSEFSS(end) S{k}.PIFSS(end) S{k}.WISFSS(end)];
        end
        R.files{w, k} = fname;
        if k == 1, caddisease1 = caddisease_INP; datatype1 = datatype_INP; dist1_1 = dist1_k; end
    end
    R.AICcs(w, :) = AICc;

    % Members must share the time grid; the pool also needs the same number of
    % draws, so a model with more draws is truncated to the smallest M.
    T = cellfun(@(c) size(c, 1), curves); M = cellfun(@(c) size(c, 2), curves);
    if any(T ~= T(1))
        error('run_ensemble:grid', 'Members have different forecast lengths (%s); same window and horizon required.', mat2str(T));
    end
    if any(M ~= M(1))
        warning('run_ensemble:draws', 'Members have different numbers of draws (%s); truncating to %d.', mat2str(M), min(M));
        for k = 1:K, curves{k} = curves{k}(:, 1:min(M)); end
    end

    % Weights: larger for the better model (see ensembleWeights).
    [weights1, winfo] = ensembleWeights(weight_type1, AICc, WISC);
    R.weights(w, :) = weights1; R.rule = winfo.rule;

    % Exact linear pool. Qens = the pool's exact quantiles at the 23 levels the scores
    % need (the weighted quantiles of ALL members' draws); forecast1 = [median LB UB]
    % for the CSV. The ensemble is SCORED from Qens (computeScoresFromQuantiles), so
    % its scores do not depend on the seed. curvesEns, a sample of trajectories from
    % the pool, is used only to draw the figure. (Formerly the ensemble was scored on
    % curvesEns, which added sampling noise of about 1% of the WIS with 6000 draws
    % per member and about 5% with 400.)
    levelsQ = [0.01 0.025 0.05:0.05:0.95 0.975 0.99];
    NensW = Nens; if isempty(NensW), NensW = size(curves{1}, 2); end
    [Qens, forecast1, ~, ~, curvesEns] = getensemble_linearpool_from_curves(weights1, curves{:}, ...
        'ProbsOut', levelsQ, 'Nens', NensW, 'Seed', seed);

    data1 = S{1}.data1(:, [1 fitIdx+1]); datalatest = S{1}.datalatest(:, [1 fitIdx+1]);   % [time, fitted variable fitIdx]
    for k = 2:K
        if ~isequal(size(S{k}.datalatest), size(S{1}.datalatest)) || ~isequal(S{k}.datalatest(:, [1 fitIdx+1]), datalatest)
            error('run_ensemble:data', 'Model %d (%s) was fitted to different data than model 1; members must share the observed series.', k, modelNames{k});
        end
    end
    nCalib = size(data1, 1);
    [perfC_e, perfF_e] = computeScoresFromQuantiles(Qens, levelsQ, datalatest, nCalib, forecastingperiod1);
    hasF = ~isempty(perfF_e);                         % false when the observations do not cover the full horizon
    R.performanceC_ens(w, :) = perfC_e;               % [MAE MSE Coverage WIS]
    if hasF
        R.performanceF_ens(w, :) = perfF_e(end, 2:5);
    end

    % Members scored through the same path (the exact quantiles of their own draws),
    % so ensemble and members are compared on the same footing. These values replace
    % the ones read from the members' files above; they agree with them up to the
    % quantile convention (the runners interpolate between draws). The weights are
    % still computed from the members' saved AICc / calibration WIS.
    for k = 1:K
        Qk = getensemble_linearpool_from_curves(1, curves{k}, 'ProbsOut', levelsQ);
        [perfC_k, perfF_k] = computeScoresFromQuantiles(Qk, levelsQ, datalatest, nCalib, forecastingperiod1);
        R.performanceC_models(w, :, k) = perfC_k;
        if ~isempty(perfF_k)
            R.performanceF_models(w, :, k) = perfF_k(end, 2:5);
        end
    end

    % --- outputs for this window ---
    suffix = strcat('-model_name-', name_ensemble, '-weight_type1-', num2str(weight_type1), '-method-', num2str(method1), ...
        '-dist-', num2str(dist1_1), '-tstart-', num2str(tstart1), '-tend-', num2str(tend1), ...
        '-calibrationperiod-', num2str(windowsize1), '-horizon-', num2str(forecastingperiod1), '-', caddisease1, '-', datatype1, '.csv');

    tv = S{1}.timevect_all; dA = S{1}.data_all(:, fitIdx); nF = length(S{1}.timevect1) + forecastingperiod1;
    if numel(tv) >= nF
        forecastdata = [tv(1:nF) dA(1:nF) forecast1];
    else
        pad = nan(nF - numel(tv), 1);
        forecastdata = [[tv; pad] [dA; pad] forecast1];
    end
    Tf = array2table(forecastdata, 'VariableNames', {'time', 'data', 'median', 'LB', 'UB'});
    writetable(Tf, fullfile(outputDir, strcat('Forecast', suffix)));

    Tc = array2table([windowsize1 R.performanceC_ens(w, :) weights1], ...
        'VariableNames', [{'Calibrationperiod', 'MAE', 'MSE', 'Coverage 95%PI', 'WIS'} strcat('w_', matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(cellstr(modelNames))))]);
    writetable(Tc, fullfile(outputDir, strcat('performance-calibration', suffix)));

    if hasF
        % Parameter-free baselines scored on the same targets with the same number
        % of draws: the forecast-hub flat baseline and a log-linear trend baseline
        % (see baselineForecast). Relative WIS = WIS / WIS_baseline; < 1 beats it.
        % NOTE: the runner's saved data1 holds the LAST BOOTSTRAP REPLICATE (it is
        % overwritten inside the bootstrap loop and never restored), so the observed
        % calibration series is taken from datalatest, which is the real data.
        calib = datalatest(1:size(data1, 1), 1:2);
        WISb = struct('flat', [], 'trend', []);
        for bt = {'flat', 'trend'}
            bdraws = baselineForecast(calib, forecastingperiod1, bt{1}, NensW, seed);
            Qb = getensemble_linearpool_from_curves(1, bdraws, 'ProbsOut', levelsQ);
            [~, perfF_b] = computeScoresFromQuantiles(Qb, levelsQ, datalatest, nCalib, forecastingperiod1);
            WISb.(bt{1}) = perfF_b(:, 5);
        end
        R.WIS_baselines(w, :) = [WISb.flat(end) WISb.trend(end)];
        R.relWIS_ens(w, :) = perfF_e(end, 5) ./ R.WIS_baselines(w, :);
        for k = 1:K
            R.relWIS_models(w, k, :) = R.performanceF_models(w, 4, k) ./ R.WIS_baselines(w, :);
        end

        % Same layout as the members' performance-forecasting file: one row per
        % horizon h, each metric averaged over horizons 1..h (as computeforecast-
        % performance / computeWIS define them), plus the baseline and relative WIS.
        perfF = [perfF_e ...
                 WISb.flat WISb.trend perfF_e(:, 5) ./ WISb.flat perfF_e(:, 5) ./ WISb.trend];
        R.performanceF_ens_byhorizon{w} = perfF;
        Tfp = array2table(perfF, 'VariableNames', {'Horizon', 'MAE', 'MSE', 'Coverage 95%PI', 'WIS', ...
            'WIS_baseline_flat', 'WIS_baseline_trend', 'relWIS_flat', 'relWIS_trend'});
        writetable(Tfp, fullfile(outputDir, strcat('performance-forecasting', suffix)));
    end

    if doPlot
        figure('Name', sprintf('%s, calibration %d', name_ensemble, windowsize1));
        for k = 1:K
            subplot(1, K + 1, k)
            plotPanel(S{k}.timevect2, curves{k}, quantile(curves{k}, [0.5 0.025 0.975], 2), S{k}.timevect1, tv, dA, forecastingperiod1, ...
                sprintf('%s (w = %.2f)', modelNames{k}, weights1(k)), caddisease1, datatype1);
        end
        subplot(1, K + 1, K + 1)
        plotPanel(S{1}.timevect2, curvesEns, forecast1, S{1}.timevect1, tv, dA, forecastingperiod1, ...
            sprintf('%s (%s)', name_ensemble, winfo.rule), caddisease1, datatype1);
    end

    fprintf('Window %d: weights [%s] (%s); ensemble calibration WIS %.3g vs members [%s]\n', windowsize1, ...
        num2str(weights1, '%.3f '), winfo.rule, R.performanceC_ens(w, 4), num2str(WISC, '%.3g '));
    if hasF
        fprintf('           horizon %d relative WIS vs flat / trend baseline: ensemble %.2f / %.2f; members [%s] / [%s]\n', ...
            forecastingperiod1, R.relWIS_ens(w, 1), R.relWIS_ens(w, 2), ...
            num2str(squeeze(R.relWIS_models(w, :, 1)), '%.2f '), num2str(squeeze(R.relWIS_models(w, :, 2)), '%.2f '));
    end
end

save(fullfile(outputDir, strcat('performancemetrics-', name_ensemble, '-weight_type1-', num2str(weight_type1), '-', caddisease1, '-', datatype1, '.mat')), 'R', '-mat');
end


function plotPanel(timevect2, draws, q, timevect1, timevect_all, data_all, horizon, ttl, caddisease1, datatype1)
% q: T x 3 [median LB UB]
plot(timevect2, draws, 'c'); hold on
plot(timevect2, q(:, 1), 'r-', 'LineWidth', 2);
plot(timevect2, q(:, 2), 'r--', 'LineWidth', 2);
plot(timevect2, q(:, 3), 'r--', 'LineWidth', 2);
plot(timevect_all, data_all, 'bo', 'LineWidth', 2);
ymax = max(q(:, 3)) * 1.5 + 1;
if horizon > 0, plot([timevect1(end) timevect1(end)], [0 ymax], 'k--', 'LineWidth', 2); end
axis([timevect1(1) timevect2(end) 0 ymax]);
xlabel('Time'); ylabel(strcat(caddisease1, {' '}, datatype1));
title(ttl); set(gca, 'FontSize', 14); set(gcf, 'color', 'white');
end
