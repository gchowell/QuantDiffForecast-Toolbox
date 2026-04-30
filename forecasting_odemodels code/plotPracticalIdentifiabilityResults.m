clear; close all; clc;

% ============================================================================
% ============================================================================
% USER SETTINGS — Change these as needed
% ============================================================================
% ============================================================================

% Options file handle (must be located in the current folder).
%   @options_forecast_PII_EXPO_r_dist1_1        (EXP, Poisson)
%   @options_forecast_PII_EXPO_r_dist1_3        (EXP, Negbin)
%   @options_forecast_PII_GGM_r_p_dist1_1       (GGM, Poisson)
%   @options_forecast_PII_GGM_r_p_dist1_3       (GGM, Negbin)
%   @options_forecast_PII_GLM_r_p_K_dist1_1     (GLM, Poisson)
%   @options_forecast_PII_GLM_r_p_K_dist1_3     (GLM, Negbin)
%   @options_forecast_PII_SIR_beta_dist1_1      (SIR S1, Poisson)
%   @options_forecast_PII_SIR_beta_dist1_3      (SIR S1, Negbin)
%   @options_forecast_PII_SIR_beta_gamma_dist1_3          (SIR S2)
%   @options_forecast_PII_SEIR_beta_dist1_3               (SEIR S1)
%   @options_forecast_PII_SEIR_beta_gamma_dist1_3         (SEIR S2)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3   (SEIR S3)
%   @options_forecast_PII_SEIRD_beta_dist1_3              (SEIRD S1)
%   @options_forecast_PII_SEIRD_beta_rho_dist1_3          (SEIRD S2)
%   @options_forecast_PII_SEIRD_beta_rho_gamma_dist1_3    (SEIRD S3)
%   @options_forecast_PII_SEIRunrep_beta_dist1_3          (SEIR_unrep S1)
%   @options_forecast_PII_SEIRunrep_beta_rho_dist1_3      (SEIR_unrep S2)
%   @options_forecast_PII_SEIRunrep_beta_rho_gamma_dist1_3 (SEIR_unrep S3)
%   @options_forecast_PII_SEIRasymp_beta0_beta1_dist1_3            (SEIR_asymp S1)
%   @options_forecast_PII_SEIRasymp_beta0_beta1_rho_dist1_3        (SEIR_asymp S2)
%   @options_forecast_PII_SEIRasymp_beta0_beta1_rho_gamma_dist1_3  (SEIR_asymp S3)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_2vars  (SEIR_vars 2vars)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_3vars  (SEIR_vars 3vars)
%   @options_forecast_PII_SEIR_beta_kappa_gamma_dist1_3_dCdt   (SEIR_vars dCdt)
options_handle = @options_forecast_PII_EXPO_r_dist1_3;

% Error structure label (for plot titles / PDF filenames):
%   'Poisson', 'Negbin5', or 'Negbin10'
error_type = 'Negbin5';

% Calibration window size — scalar or array
%   scalar:  windowsize1 = 20;
%   array:   windowsize1 = 20:10:100;
windowsize1 = 30:10:50;

% Number of replicates
num_replicates = 20;

% Whether to run the computation (true) or skip to plotting (false)
run_flag = true;

% Plot type (after computation or from existing output):
%   'PII'     — PII median + 95% CI vs calibration period
%   'CI_grid' — 95% CI caterpillar plot per parameter
%   'both'    — generate both PII and CI_grid plots
%   'none'    — no plotting (run only)
plot_type = 'both';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================================================================
% ============================================================================
% END USER SETTINGS — Change these as needed
% ============================================================================
% ============================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================================================================
% EXTRACT MODEL CONFIGURATION FROM OPTIONS FILE
% ============================================================================

[~, ~, ~, ~, ~, ~, model_INP, params_INP, ~, ~, ~, ~, ~, ~, ~] = options_handle();

model_name     = model_INP.name;   % full name for file loading
parts = strsplit(model_name, '-PII');
model_display  = parts{1};         % short name for titles (EXP, GGM, SIR, SEIR, ...)

params_label   = params_INP.label;
params_fixed   = params_INP.fixed;
params_initial = params_INP.initial;

estimated_indices = find(params_fixed == 0);
num_estimated     = length(estimated_indices);
true_values       = params_initial(estimated_indices);

fprintf('Model      : %s\n', model_display);
fprintf('Full name  : %s\n', model_name);
fprintf('Error type : %s\n', error_type);
fprintf('Window sizes: %s\n', mat2str(windowsize1));
fprintf('Estimated parameters: ');
for k = 1:num_estimated
    fprintf('%s ', params_label{estimated_indices(k)});
end
fprintf('\n\n');

basePath = fileparts(mfilename('fullpath'));
output_root_Dir = fullfile(basePath, 'output');
analysis_folder = make_analysis_folder_name(model_name, model_display, error_type, windowsize1);
output_Dir = fullfile(output_root_Dir, analysis_folder);

if ~exist(output_Dir, 'dir')
    mkdir(output_Dir);
end

fprintf('Output folder: %s\n\n', output_Dir);

% ============================================================================
% RUN COMPUTATION (parallel across replicates)
% ============================================================================

if run_flag
    num_ws = length(windowsize1);

    % Build list of all (windowsize, replicate) tasks
    total_tasks = num_ws * num_replicates;
    task_ws  = zeros(1, total_tasks);
    task_rep = zeros(1, total_tasks);
    idx = 0;
    for w = 1:num_ws
        for r = 1:num_replicates
            idx = idx + 1;
            task_ws(idx)  = windowsize1(w);
            task_rep(idx) = r;
        end
    end

    % Filter out tasks whose output already exists
    run_mask = true(1, total_tasks);
    for t = 1:total_tasks

        outfile = result_file_path(output_Dir, task_rep(t), task_ws(t), model_name, false);

        if isfile(outfile)
            fprintf('Output exists: %s - skipping.\n', outfile);
            run_mask(t) = false;
        end

    end
    run_idx = find(run_mask);
    num_to_run = length(run_idx);

    if num_to_run > 0
        fprintf('\nRunning %d tasks (parfor) ...\n', num_to_run);

        % Extract the tasks to run
        run_ws  = task_ws(run_idx);
        run_rep = task_rep(run_idx);

        parfor t = 1:num_to_run
            rng(run_rep(t));
            fprintf('Running: replicate %d, windowsize %d ...\n', ...
                run_rep(t), run_ws(t));
            Run_PracticalIndentifiability_ODEModel(...
                options_handle, run_ws(t), 5, 1, run_rep(t), output_Dir);
        end
    end

    fprintf('\nComputation complete.\n\n');
end

% ============================================================================
% PLOTTING (reads from ./output/)
% ============================================================================

%output_dir = fullfile('.', 'output');

close all;

if ismember(plot_type, {'PII', 'both'})
    plot_PII_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, windowsize1, num_replicates);
end

if ismember(plot_type, {'CI_grid', 'both'})
    plot_CI_grid_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, true_values, ...
        windowsize1, num_replicates);
end

if strcmp(plot_type, 'none')
    fprintf('Plotting skipped.\n');
end

% ============================================================================
% FUNCTION: plot_PII_figure
% ============================================================================
function plot_PII_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, windowsize1, num_replicates)

    num_estimated = length(estimated_indices);
    num_ws = length(windowsize1);

    fs_ylabel = 18;
    fs_title  = 18;
    fs_legend = 18;
    fs_tick   = 18;

    colors_list = {'r', 'b', [0 0.6 0], [0.8 0.4 0], [0.5 0 0.5]};

    if num_ws == 1
        % --- Single window size: boxplot per parameter ---
        figure('Position', [100, 100, 350*num_estimated, 300]);
        tiledlayout(1, num_estimated, 'Padding', 'compact', 'TileSpacing', 'compact');

        for p = 1:num_estimated
            clr = colors_list{min(p, length(colors_list))};
            SCI_vals = NaN(num_replicates, 1);
            count = 0;
            for r = 1:num_replicates
                fname = result_file_path(output_Dir, r, windowsize1, model_name, true);
                if isfile(fname)
                    S = load(fname, 'SCI');
                    SCI_vals(r) = S.SCI(1, estimated_indices(p));
                    count = count + 1;
                end
            end
            valid = SCI_vals(~isnan(SCI_vals));
            if ~isempty(valid)
                lb = quantile(valid, 0.025); med = quantile(valid, 0.5); ub = quantile(valid, 0.975);
            else
                lb = NaN; med = NaN; ub = NaN;
            end
            fprintf('%s | param=%s | T=%d : %d files, LB=%.6f, Median=%.6f, UB=%.6f\n', ...
                error_type, params_label{estimated_indices(p)}, windowsize1, count, lb, med, ub);

            nexttile;
            boxplot(valid, 'Labels', {sprintf('T=%d', windowsize1)});
            hold on;
            if ~isempty(valid)
                plot(1, med, 'o', 'Color', clr, 'MarkerSize', 8, ...
                    'MarkerFaceColor', clr, 'LineWidth', 2);
            end
            if ub > 1
                yline(1, 'k--', 'LineWidth', 1.2);
            end
            ylabel(sprintf('PII (%s)', params_label{estimated_indices(p)}), ...
                'Interpreter', 'latex', 'FontSize', fs_ylabel);
            if p == 1
                title(sprintf('%s, %s', model_display, error_type), ...
                    'Interpreter', 'latex', 'FontSize', fs_title);
            end
            set(gca, 'FontSize', fs_tick, 'TickLabelInterpreter', 'latex');
            grid off; box on; ylim([0 inf]);
        end

        outname = 'plot_box_param';

    else
        % --- Multiple window sizes: line plot (median + 95% CI) ---
        figure('Position', [100, 100, 400, 350*num_estimated]);
        tiledlayout(num_estimated, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

        for p = 1:num_estimated
            clr = colors_list{min(p, length(colors_list))};
            SCIss = NaN(num_replicates, num_ws);
            for w = 1:num_ws
                ws = windowsize1(w);
                count = 0;
                for r = 1:num_replicates
                    fname = result_file_path(output_Dir, r, ws, model_name, true);
                    if isfile(fname)
                        S = load(fname, 'SCI');
                        SCIss(r, w) = S.SCI(1, estimated_indices(p));
                        count = count + 1;
                    end
                end
                valid = SCIss(~isnan(SCIss(:,w)), w);
                if ~isempty(valid)
                    lb = quantile(valid, 0.025); med = quantile(valid, 0.5); ub = quantile(valid, 0.975);
                else
                    lb = NaN; med = NaN; ub = NaN;
                end
                fprintf('%s | param=%s | T=%d : %d files, LB=%.6f, Median=%.6f, UB=%.6f\n', ...
                    error_type, params_label{estimated_indices(p)}, ws, count, lb, med, ub);
            end

            ax = nexttile;
            med_vals = max(quantile(SCIss, 0.5, 1), 0);
            lb_vals  = max(quantile(SCIss, 0.025, 1), 0);
            ub_vals  = max(quantile(SCIss, 0.975, 1), 0);

            plot(windowsize1, med_vals, '-o', 'Color', clr, ...
                'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', clr);
            hold on;
            plot(windowsize1, lb_vals, 'k--', 'LineWidth', 1.5);
            plot(windowsize1, ub_vals, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            if any(ub_vals > 1)
                yline(1, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            end
            xlim([windowsize1(1) windowsize1(end)]);
            xticks(windowsize1);
            ylabel(sprintf('PII (%s)', params_label{estimated_indices(p)}), ...
                'Interpreter', 'latex', 'FontSize', fs_ylabel);
            if p == 1
                title(sprintf('%s, %s', model_display, error_type), ...
                    'Interpreter', 'latex', 'FontSize', fs_title);
            end
            legend({'Median', '95\% CI'}, ...
                'Interpreter', 'latex', 'Location', 'best', 'FontSize', fs_legend);
            set(gca, 'FontSize', fs_tick, 'TickLabelInterpreter', 'latex');
            grid off; box on; ylim([0 inf]);
            axs(p) = ax; %#ok<AGROW>
        end
        if num_estimated > 1
            linkaxes(axs, 'x');
        end

        outname = 'PII_param';
    end

    set(gcf, 'Color', 'white', 'Renderer', 'painters');
    set(gcf, 'Units', 'centimeters');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [pos(3) pos(4)], ...
        'PaperPosition', [0 0 pos(3) pos(4)]);
    outfile = fullfile(output_Dir, outname);
    print(gcf, outfile, '-dpdf');
    fprintf('Saved: %s.pdf\n', outfile);
end

% ============================================================================
% FUNCTION: plot_CI_grid_figure
%   Rows = calibration periods (windowsize1 values)
%   Each cell: CI caterpillar plot — gray if CI contains true value, red otherwise
%   Blue dashed vertical line at true parameter value
%   One figure per estimated parameter
% ============================================================================
function plot_CI_grid_figure(output_Dir, error_type, model_name, model_display, ...
        params_label, estimated_indices, true_values, ...
        windowsize1, num_replicates)

    num_estimated = length(estimated_indices);
    num_ws = length(windowsize1);

    fs_title = 14;
    fs_label = 14;
    fs_tick  = 12;

    for param_i = 1:num_estimated
        param_idx = estimated_indices(param_i);

        fig = figure('Position', [50, 50, 350, 250*num_ws]);
        tiledlayout(num_ws, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

        for w = 1:num_ws
            ws = windowsize1(w);

            estimates = [];
            LBs = [];
            UBs = [];

            for r = 1:num_replicates
                fname = result_file_path(output_Dir, r, ws, model_name, true);
                if isfile(fname)
                    S = load(fname, 'paramss');
                    col_median = 1 + (param_idx-1)*3;
                    col_LB     = 2 + (param_idx-1)*3;
                    col_UB     = 3 + (param_idx-1)*3;
                    estimates(end+1) = S.paramss(1, col_median); %#ok<AGROW>
                    LBs(end+1)       = S.paramss(1, col_LB);     %#ok<AGROW>
                    UBs(end+1)       = S.paramss(1, col_UB);     %#ok<AGROW>
                end
            end

            nexttile;
            hold on;

            num_sims = length(estimates);
            for sim_i = 1:num_sims
                contains_true = (LBs(sim_i) <= true_values(param_i)) && ...
                                (true_values(param_i) <= UBs(sim_i));
                if contains_true
                    plot([LBs(sim_i) UBs(sim_i)], [sim_i sim_i], '-', ...
                        'Color', [0.6 0.6 0.6], 'LineWidth', 1);
                else
                    plot([LBs(sim_i) UBs(sim_i)], [sim_i sim_i], 'r-', ...
                        'LineWidth', 1);
                end
                plot(estimates(sim_i), sim_i, 'k.', 'MarkerSize', 8);
            end

            yL = [0 max(num_sims,1)+1];
            plot([true_values(param_i) true_values(param_i)], yL, 'b--', 'LineWidth', 2);
            ylim(yL);
            grid off;
            set(gca, 'FontSize', fs_tick, 'TickLabelInterpreter', 'latex');

            if w == 1
                title(sprintf('%s, %s, %s', model_display, error_type, ...
                    params_label{param_idx}), 'Interpreter', 'latex', 'FontSize', fs_title);
            end
            ylabel(sprintf('$T = %d$', ws), ...
                'Interpreter', 'latex', 'FontSize', fs_label);
            if w == num_ws
                xlabel('Estimate', 'Interpreter', 'latex', 'FontSize', fs_label);
            end
        end

        set(fig, 'Color', 'white', 'Renderer', 'painters');
        set(fig, 'Units', 'centimeters');
        pos = get(fig, 'Position');
        set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [pos(3) pos(4)], ...
            'PaperPosition', [0 0 pos(3) pos(4)]);
        outname = sprintf('CI_grid_param_%d', param_idx);
        outfile = fullfile(output_Dir, outname);
        print(fig, outfile, '-dpdf');
        fprintf('Saved: %s.pdf\n', outfile);
    end
end

function fname = result_file_path(output_Dir, replicate_id, windowsize1, model_name, use_legacy)
    fname = fullfile(output_Dir, ...
        sprintf('results_%s_rep%d.mat', window_file_token(windowsize1), replicate_id));

    if use_legacy && ~isfile(fname)
        legacy_Dir = fileparts(output_Dir);
        legacy_name = sprintf('results-replicate-%d-model_name-%s-calibrationperiod-%d.mat', ...
            replicate_id, model_name, windowsize1);
        legacy_fname = fullfile(legacy_Dir, legacy_name);
        if isfile(legacy_fname)
            fname = legacy_fname;
        end
    end
end

function folder_name = make_analysis_folder_name(model_name, model_display, error_type, windowsize1)
    scenario = scenario_from_model_name(model_name);
    window_label = window_folder_label(windowsize1);

    if isempty(scenario) || strcmp(scenario, 'r')
        folder_name = sprintf('%s_%s_%s', model_display, error_type, window_label);
    else
        folder_name = sprintf('%s_%s_%s_%s', model_display, scenario, error_type, window_label);
    end

    folder_name = regexprep(folder_name, '[<>:"/\\|?*]', '_');
    folder_name = regexprep(folder_name, '\s+', '_');
    folder_name = regexprep(folder_name, '_+', '_');
end

function scenario = scenario_from_model_name(model_name)
    parts = strsplit(model_name, '-PII-');
    if numel(parts) < 2
        scenario = '';
        return
    end

    scenario = regexprep(parts{2}, '-dist1-\d+', '');
    scenario = regexprep(scenario, '^-|-$', '');
end

function label = window_folder_label(windowsize1)
    if numel(windowsize1) == 1
        label = sprintf('W%s', number_file_token(windowsize1));
    else
        diffs = diff(windowsize1);
        if all(diffs == diffs(1))
            label = sprintf('W%s-%s-%s', number_file_token(windowsize1(1)), ...
                number_file_token(diffs(1)), number_file_token(windowsize1(end)));
        else
            pieces = arrayfun(@number_file_token, windowsize1, 'UniformOutput', false);
            label = ['W' strjoin(pieces, '-')];
        end
    end
end

function token = window_file_token(windowsize1)
    token = sprintf('W%s', number_file_token(windowsize1));
end

function token = number_file_token(value1)
    token = strrep(sprintf('%g', value1), '.', 'p');
end
