function make_readme_figures(rerun)
%MAKE_README_FIGURES  Generate the additional example figures shown in README.md.
%
%   make_readme_figures          run what is missing, then draw and save the figures
%   make_readme_figures(true)    discard the cached results and rerun everything
%
%   Place this file in 'forecasting_odemodels code' and run it with that folder
%   on the MATLAB path. It writes four PNG files to <repository>/docs/images:
%
%     forecast_by_calibration_length.png  10-day SEIR forecasts made from the
%                                         first 14, 17 and 20 observations
%     noise_model_comparison.png          the same forecast under least squares
%                                         with normal noise, Poisson likelihood
%                                         and negative-binomial likelihood
%     baseline_comparison.png             SEIR forecast against the flat and
%                                         trend baselines, with WIS by horizon
%     profile_likelihood_beta.png         profile likelihood of beta
%
%   The script calls the shipped options files unchanged (five forecast runs and
%   one profile), so the runs use the bootstrap settings in those files. Results
%   are summarized into output/readme_figures_cache.mat after every run; a
%   second call only redraws the figures. The runners call "close all", so all
%   figures are drawn after the last run.
%
%   Data: daily influenza counts, San Francisco 1918 (input/curve-flu1918SF.txt).

if nargin < 1 || isempty(rerun), rerun = false; end

% <============================================================================>
% <=========================== Settings =======================================>
% <============================================================================>

H    = 10;      % forecast horizon (days)
seed = 1;       % global RNG seed set before every run
W0   = 17;      % calibration length of the quick-start example

% tag, forecast options, calibration length
runs = { ...
    'nb14',   @options_forecast_SEIR_flu1918_dist1_3, 14; ...
    'nb17',   @options_forecast_SEIR_flu1918_dist1_3, W0; ...
    'nb20',   @options_forecast_SEIR_flu1918_dist1_3, 20; ...
    'lsq17',  @options_forecast_SEIR_flu1918_dist1_0, W0; ...
    'pois17', @options_forecast_SEIR_flu1918_dist1_1, W0};

blue   = [0.18 0.43 0.70];
orange = [0.80 0.36 0.04];
gray   = [0.45 0.47 0.50];
dark   = [0.15 0.15 0.15];

% <============================================================================>
% <=========================== Folders and data ===============================>
% <============================================================================>

codeDir = fileparts(which('Run_Forecasting_ODEModel'));
if isempty(codeDir)
    error('make_readme_figures:path', ...
        'Run_Forecasting_ODEModel is not on the MATLAB path. Add ''forecasting_odemodels code'' to the path first.');
end
imgDir = fullfile(fileparts(codeDir), 'docs', 'images');
outDir = fullfile(codeDir, 'output');
if ~exist(imgDir, 'dir'), mkdir(imgDir); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
cacheFile = fullfile(outDir, 'readme_figures_cache.mat');

data = load(fullfile(codeDir, 'input', 'curve-flu1918SF.txt'));

% <============================================================================>
% <=========================== Runs (cached) ==================================>
% <============================================================================>

R = struct();
PLres = [];
if ~rerun && exist(cacheFile, 'file') == 2
    S = load(cacheFile);
    if isfield(S, 'R'),     R = S.R;         end
    if isfield(S, 'PLres'), PLres = S.PLres; end
end

for k = 1:size(runs, 1)
    tag = runs{k, 1};
    if isfield(R, tag), continue; end
    fprintf('\n===== make_readme_figures: forecast run %d of %d (%s) =====\n', k, size(runs, 1), tag);
    rng(seed, 'twister');
    [~, ~, ~, ~, draws, calib] = Run_Forecasting_ODEModel(runs{k, 2}, 1, 1, runs{k, 3}, H);
    R.(tag) = summarizeDraws(draws, calib, data, H);
    save(cacheFile, 'R', 'PLres');
end

if isempty(PLres)
    fprintf('\n===== make_readme_figures: profile likelihood of beta =====\n');
    rng(seed, 'twister');
    plOpts = struct('paramIndex', 1, 'alpha', 0.05, 'numstartpoints', 3, 'printscreen', 0);
    PL = Run_ProfileLikelihood_ODEModel(@options_fit_SEIR_flu1918_dist1_3, 1, W0, plOpts);
    p = PL.param(1);
    PLres = struct('theta', p.theta, 'D', p.D, 'mle', p.mle, 'CI', p.CI, ...
        'threshold', PL.threshold, 'alpha', PL.alpha, 'status', p.status);
    save(cacheFile, 'R', 'PLres');
end

% Baselines on the quick-start targets: no fitting, deterministic for a given seed.
calib0 = data(1:W0, :);
Bflat  = summarizeDraws(baselineForecast(calib0, H, 'flat',  6000, 1), calib0, data, H);
Btrend = summarizeDraws(baselineForecast(calib0, H, 'trend', 6000, 1), calib0, data, H);

close all

% <============================================================================>
% <============== Figure 1: forecasts as calibration data accrue ==============>
% <============================================================================>

tags = {'nb14', 'nb17', 'nb20'};
lastRow = max(cellfun(@(t) numel(R.(t).time), tags));
ytop = commonTop(R, tags, max(data(1:lastRow, 2)));

fig = figure('Color', 'white', 'Position', [60 60 1400 430]);
tl = tiledlayout(fig, 1, numel(tags), 'Padding', 'compact', 'TileSpacing', 'compact');
for k = 1:numel(tags)
    s = R.(tags{k});
    ax = nexttile(tl);
    h = forecastPanel(ax, s, blue, gray, dark);
    xlim(ax, [data(1, 1) - 0.5, data(lastRow, 1) + 0.5]);
    ylim(ax, [0 ytop]);
    title(ax, {sprintf('Calibration: days %g-%g', s.time(1), s.time(s.W)), scoreLine(s)}, 'FontWeight', 'normal');
    xlabel(ax, 'Time (days)');
    if k == 1
        ylabel(ax, 'New cases per day');
        legend(ax, h, {'95% PI', '50% PI', 'Predictive median', 'Calibration data', 'Held-out data'}, ...
            'Location', 'northwest', 'Box', 'off');
    end
end
title(tl, sprintf('SEIR model, negative-binomial likelihood: %d-day forecasts from three calibration lengths', H), 'FontSize', 13);
saveFigure(fig, fullfile(imgDir, 'forecast_by_calibration_length.png'));

% <============================================================================>
% <============== Figure 2: estimation method and noise model =================>
% <============================================================================>

tags   = {'lsq17', 'pois17', 'nb17'};
labels = {'Least squares, normal noise (method 0, dist 0)', ...
          'Poisson likelihood (method 1)', ...
          'Negative-binomial likelihood (method 3)'};
lastRow = W0 + H;
ytop = commonTop(R, tags, max(data(1:lastRow, 2)));

fig = figure('Color', 'white', 'Position', [60 60 1400 430]);
tl = tiledlayout(fig, 1, numel(tags), 'Padding', 'compact', 'TileSpacing', 'compact');
for k = 1:numel(tags)
    s = R.(tags{k});
    ax = nexttile(tl);
    h = forecastPanel(ax, s, blue, gray, dark);
    xlim(ax, [data(1, 1) - 0.5, data(lastRow, 1) + 0.5]);
    ylim(ax, [0 ytop]);
    title(ax, {labels{k}, scoreLine(s)}, 'FontWeight', 'normal');
    xlabel(ax, 'Time (days)');
    if k == 1
        ylabel(ax, 'New cases per day');
        legend(ax, h, {'95% PI', '50% PI', 'Predictive median', 'Calibration data', 'Held-out data'}, ...
            'Location', 'northwest', 'Box', 'off');
    end
end
title(tl, sprintf('Same SEIR model and data (first %d observations), three observation models', W0), 'FontSize', 13);
saveFigure(fig, fullfile(imgDir, 'noise_model_comparison.png'));

% <============================================================================>
% <============== Figure 3: model against the baselines =======================>
% <============================================================================>

M = R.nb17;
f = (W0 + 1):(W0 + H);                 % held-out rows
members = {M, Btrend, Bflat};
cols    = {blue, orange, gray};
names   = {'SEIR model', 'Trend baseline', 'Flat baseline'};

fig = figure('Color', 'white', 'Position', [60 60 1100 430]);
tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax = nexttile(tl); hold(ax, 'on');
hm = gobjects(1, numel(members));
for m = 1:numel(members)
    s = members{m};
    t = s.time(f);
    fill(ax, [t; flipud(t)], [max(s.lo95(f), 1); flipud(max(s.hi95(f), 1))], cols{m}, ...
        'FaceAlpha', 0.14, 'EdgeColor', 'none');
end
for m = 1:numel(members)
    s = members{m};
    hm(m) = plot(ax, s.time(f), max(s.med(f), 1), '-', 'Color', cols{m}, 'LineWidth', 2);
end
c = 1:W0;
hc = plot(ax, M.time(c), max(M.obs(c), 1), 'o', 'MarkerSize', 6, 'MarkerFaceColor', dark, 'MarkerEdgeColor', 'w');
hf = plot(ax, M.time(f), max(M.obs(f), 1), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', dark, 'LineWidth', 1.2);
xline(ax, M.time(W0) + 0.5, '--', 'Color', gray);
set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
xlim(ax, [M.time(1) - 0.5, M.time(end) + 0.5]);
xlabel(ax, 'Time (days)');
ylabel(ax, 'New cases per day (log scale)');
title(ax, 'Predictive medians and 95% prediction intervals', 'FontWeight', 'normal');
legend(ax, [hm hc hf], [names, {'Calibration data', 'Held-out data'}], 'Location', 'northwest', 'Box', 'off');

ax = nexttile(tl); hold(ax, 'on');
hw = gobjects(1, numel(members));
for m = 1:numel(members)
    hw(m) = plot(ax, (1:H)', members{m}.wis, '-o', 'Color', cols{m}, 'LineWidth', 2, ...
        'MarkerSize', 5, 'MarkerFaceColor', cols{m});
end
set(ax, 'YScale', 'log', 'FontSize', 11, 'Box', 'on');
xlim(ax, [0.5, H + 0.5]);
xlabel(ax, 'Forecast horizon h (days)');
ylabel(ax, 'Mean WIS over leads 1 to h (log scale)');
relTrend = M.wis(end) / Btrend.wis(end);
relFlat  = M.wis(end) / Bflat.wis(end);
title(ax, sprintf('Relative WIS at h = %d: %.2f vs trend, %.2f vs flat', H, relTrend, relFlat), 'FontWeight', 'normal');
legend(ax, hw, names, 'Location', 'northwest', 'Box', 'off');

title(tl, sprintf('SEIR forecast against parameter-free baselines (first %d observations)', W0), 'FontSize', 13);
saveFigure(fig, fullfile(imgDir, 'baseline_comparison.png'));

% <============================================================================>
% <============== Figure 4: profile likelihood of beta ========================>
% <============================================================================>

ok = isfinite(PLres.theta) & isfinite(PLres.D);
[th, order] = sort(PLres.theta(ok));
D = PLres.D(ok);
D = D(order);
ytop = max(1.5 * PLres.threshold, min(max([D; 0]) * 1.1, 6 * PLres.threshold));

fig = figure('Color', 'white', 'Position', [60 60 660 430]);
ax = axes(fig); hold(ax, 'on');
if all(isfinite(PLres.CI))
    fill(ax, PLres.CI([1 2 2 1]), [0 0 ytop ytop], blue, 'FaceAlpha', 0.10, 'EdgeColor', 'none');
end
plot(ax, th, D, '-o', 'Color', blue, 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', blue);
yline(ax, PLres.threshold, '--', 'Color', gray, 'LineWidth', 1);
plot(ax, PLres.mle, 0, 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
if numel(th) > 1, xlim(ax, [th(1) th(end)]); end
ylim(ax, [0 ytop]);
set(ax, 'FontSize', 11, 'Box', 'on');
xlabel(ax, '\beta', 'Interpreter', 'tex');
ylabel(ax, '-2 \Delta log L', 'Interpreter', 'tex');
level = 100 * (1 - PLres.alpha);
title(ax, {sprintf('Profile likelihood of \\beta, first %d observations', W0), ...
    sprintf('MLE %.3g, %.0f%% interval [%.3g, %.3g]', PLres.mle, level, PLres.CI(1), PLres.CI(2))}, ...
    'FontWeight', 'normal', 'Interpreter', 'tex');
text(ax, 0.98, 0.95, sprintf('dashed line: \\chi^2_1(%.2f) = %.2f', 1 - PLres.alpha, PLres.threshold), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'Interpreter', 'tex');
saveFigure(fig, fullfile(imgDir, 'profile_likelihood_beta.png'));

% <============================================================================>
% <=========================== Console summary ================================>
% <============================================================================>

fprintf('\n%-10s %-12s %-22s %s\n', 'run', 'calibration', '95% PI coverage (%)', sprintf('WIS, leads 1-%d', H));
allTags = runs(:, 1)';
for k = 1:numel(allTags)
    s = R.(allTags{k});
    fprintf('%-10s %-12d %-22.0f %.2f\n', allTags{k}, s.W, s.coverage, s.wis(end));
end
fprintf('%-10s %-12d %-22.0f %.2f\n', 'flat',  Bflat.W,  Bflat.coverage,  Bflat.wis(end));
fprintf('%-10s %-12d %-22.0f %.2f\n', 'trend', Btrend.W, Btrend.coverage, Btrend.wis(end));
fprintf('Profile of beta: MLE %.4g, %.0f%% interval [%.4g, %.4g], status: %s\n', ...
    PLres.mle, level, PLres.CI(1), PLres.CI(2), PLres.status);
fprintf('\nFigures saved in %s\n', imgDir);

end

% <============================================================================>
% <=========================== Local functions ================================>
% <============================================================================>

function s = summarizeDraws(draws, calib, data, H)
% Quantiles, WIS by horizon and 95% PI coverage of a (W+H) x Ndraws matrix of
% predictive draws for the first observed series. Quantiles are clipped at zero,
% as in the runner's exports.
W = size(calib, 1);
rows = 1:(W + H);
if size(draws, 1) < W + H
    error('make_readme_figures:draws', 'Expected %d rows of predictive draws, found %d.', W + H, size(draws, 1));
end
if size(data, 1) < W + H
    error('make_readme_figures:data', 'The data file has %d rows; %d are needed to score this forecast.', size(data, 1), W + H);
end
draws = draws(rows, :);
q = max(quantile(draws, [0.025 0.25 0.5 0.75 0.975], 2), 0);

s.W    = W;
s.time = data(rows, 1);
s.obs  = data(rows, 2);
s.lo95 = q(:, 1);
s.lo50 = q(:, 2);
s.med  = q(:, 3);
s.hi50 = q(:, 4);
s.hi95 = q(:, 5);

[~, wis] = computeWIS(calib(:, 1:2), data(:, 1:2), draws, H);
s.wis = wis(:, 2);                       % mean WIS over leads 1..h, h = 1..H

f = (W + 1):(W + H);
s.coverage = 100 * mean(s.obs(f) >= s.lo95(f) & s.obs(f) <= s.hi95(f));
end

function h = forecastPanel(ax, s, col, gray, dark)
% Fit and forecast of one run: prediction bands, median, calibration and held-out data.
hold(ax, 'on');
t = s.time;
c = 1:s.W;
f = (s.W + 1):numel(t);
h = gobjects(1, 5);
h(1) = fill(ax, [t; flipud(t)], [s.lo95; flipud(s.hi95)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
h(2) = fill(ax, [t; flipud(t)], [s.lo50; flipud(s.hi50)], col, 'FaceAlpha', 0.30, 'EdgeColor', 'none');
h(3) = plot(ax, t, s.med, '-', 'Color', col, 'LineWidth', 2);
h(4) = plot(ax, t(c), s.obs(c), 'o', 'MarkerSize', 6, 'MarkerFaceColor', dark, 'MarkerEdgeColor', 'w');
h(5) = plot(ax, t(f), s.obs(f), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', dark, 'LineWidth', 1.2);
xline(ax, t(s.W) + 0.5, '--', 'Color', gray);
set(ax, 'FontSize', 11, 'Box', 'on');
end

function txt = scoreLine(s)
txt = sprintf('95%% PI coverage %.0f%%, WIS %.1f (held-out days)', s.coverage, s.wis(end));
end

function ytop = commonTop(R, tags, maxObs)
% Shared upper y-limit: the largest 95% bound across panels, capped at 2.5 times
% the largest plotted observation so one wide interval does not flatten the rest.
top = 0;
for k = 1:numel(tags)
    top = max(top, max(R.(tags{k}).hi95));
end
ytop = 1.1 * max(maxObs, min(top, 2.5 * maxObs));
end

function saveFigure(fig, file)
drawnow;
try
    exportgraphics(fig, file, 'Resolution', 200);
catch
    print(fig, file, '-dpng', '-r200');
end
fprintf('Saved: %s\n', file);
end
