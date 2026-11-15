function results = qdfSmokeTest2(tests, cleanup)
%QDFSMOKETEST2  Second set of MATLAB checks for QuantDiffForecast: the paths qdfSmokeTest does not cover.
%
%   qdfSmokeTest2                 run all ten tests, then delete the test files
%   qdfSmokeTest2([1 2 3])        run only tests 1 to 3
%   qdfSmokeTest2(1:10, false)    run all tests and keep the test files and outputs
%   results = qdfSmokeTest2(...)  also return the summary as a cell array
%
%   Run it from the toolbox code folder. It needs qdfTestOptions.m in the same
%   folder. Each test uses a quick copy of a shipped options file (20 bootstrap
%   fits, its own model name), so your real results in 'output' are not touched.
%
%   Tests (1918 influenza data, window of 17 points unless stated)
%      1  Poisson likelihood (method 1)
%      2  negative binomial, variance = mean + alpha*mean^2 (method 4)
%      3  negative binomial, variance = mean + alpha*mean^d (method 5)
%      4  least squares with normal bootstrap noise (method 0, dist1 0)
%      5  least squares with Poisson bootstrap noise (method 0, dist1 1)
%      6  Poisson fallback: least squares with NB noise (dist1 2) on a smooth
%         series whose dispersion is below 1; file names must stay '-dist-2-'
%         and plotFit_ODEModel must find them
%      7  initial value estimated (params.fixI0 = 0)
%      8  two fitted series (an increment series and a level series)
%      9  rolling windows: three window starts, then plotForecast_ODEModel
%     10  identifiability study (2 replicates) and summarizeIdentifiability
%
%   For tests 1-5, 7 and 8 the script checks the number of estimated quantities
%   and the AICc penalty term, which depends on the effective sample size, and
%   the noise model recorded in noiseLog.
%
%   The console output of the runners is hidden; only a one-line verdict per
%   test and a summary are printed. A FAIL line carries the error message or the
%   value that was not as expected: send those lines back for diagnosis.

if nargin < 1 || isempty(tests),   tests = 1:10;    end
if nargin < 2 || isempty(cleanup), cleanup = true;  end

if exist('qdfTestOptions', 'file') ~= 2
    error('qdfSmokeTest2:helper', 'qdfTestOptions.m must be in the current folder or on the path.');
end
if exist('Run_Fit_ODEModel', 'file') ~= 2
    error('qdfSmokeTest2:folder', 'Run this from the toolbox code folder (Run_Fit_ODEModel.m not found).');
end
if ~isnumeric(tests) || any(tests(:) < 1) || any(tests(:) > 10) || any(tests(:) ~= fix(tests(:)))
    error('qdfSmokeTest2:tests', 'tests must be whole numbers between 1 and 10.');
end

names = { ...
    'Poisson likelihood (method 1)', ...
    'negative binomial, mean^2 variance (method 4)', ...
    'negative binomial, mean^d variance (method 5)', ...
    'least squares + normal noise (dist1 0)', ...
    'least squares + Poisson noise (dist1 1)', ...
    'Poisson fallback keeps -dist-2- file names', ...
    'initial value estimated (fixI0 = 0)', ...
    'two fitted series', ...
    'rolling windows (3 starts) + plotForecast', ...
    'identifiability study + summary'};
results = cell(0, 3);
fit3 = 'options_fit_SEIR_flu1918_dist1_3';        % negative binomial example
fit0 = 'options_fit_SEIR_flu1918_dist1_0';        % least-squares example
fc3  = 'options_forecast_SEIR_flu1918_dist1_3';
pii3 = 'options_forecast_PII_EXPO_r_dist1_3';
mLine = 'method1 = 3; % Parameter estimation method';
dLine = 'dist1 = 0; % Error structure type';
outPII = fullfile(pwd, 'output', 'qdf_test_PII');

for t = tests(:)'
    fprintf('Test %d: %s ... ', t, names{t});
    t0 = tic;
    try
        switch t
            case {1, 2, 3}
                meth = [1 4 5];  meth = meth(t);
                kExp = [1 2 3];  kExp = kExp(t);                 % beta (+ alpha) (+ d)
                opt = sprintf('options_test_u%d', t);  mdl = sprintf('SEIRu%d', t);
                qdfQuiet(@() qdfTestOptions(fit3, opt, 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', ['model.name = ''' mdl ''';'], ...
                    mLine, sprintf('method1 = %d; %% Parameter estimation method', meth)));
                rng(1, 'twister');
                [~, A] = evalc(sprintf('Run_Fit_ODEModel(@%s, 1, 1, 17)', opt));
                L = qdfNoiseLog(['bootstraps-ODEModel-*' mdl '*.mat']);
                pen = qdfPenalty(kExp, 16);                       % 17 points, first one fixed
                ok = isreal(A) && isfinite(A(1, 2)) && A(1, 5) == kExp && abs(A(1, 4) - pen) < 1e-6 && ...
                     L(1, 2) == meth && L(1, 3) == meth && L(1, 4) > 0;
                if t == 1, ok = ok && L(1, 4) == 1; end
                if t == 3, ok = ok && L(1, 5) >= 0.6 && L(1, 5) <= 5; end
                detail = sprintf('AICc = %.4g, numparams = %d (expected %d), penalty = %.4f (expected %.4f); noiseLog = [%s]', ...
                    A(1, 2), A(1, 5), kExp, A(1, 4), pen, num2str(L(1, :), '%.4g '));

            case {4, 5}
                dval = t - 4;                                     % 0 = normal, 1 = Poisson
                opt = sprintf('options_test_u%d', t);  mdl = sprintf('SEIRu%d', t);
                if dval == 0
                    qdfQuiet(@() qdfTestOptions(fit0, opt, 'B = 300;', 'B = 20;', ...
                        'model.name = ''SEIR model'';', ['model.name = ''' mdl ''';']));
                else
                    qdfQuiet(@() qdfTestOptions(fit0, opt, 'B = 300;', 'B = 20;', ...
                        'model.name = ''SEIR model'';', ['model.name = ''' mdl ''';'], ...
                        dLine, 'dist1 = 1; % Error structure type'));
                end
                rng(1, 'twister');
                [~, A] = evalc(sprintf('Run_Fit_ODEModel(@%s, 1, 1, 17)', opt));
                L = qdfNoiseLog(['bootstraps-ODEModel-*' mdl '*.mat']);
                part2 = 16 * log(2 * pi) + 16 + qdfPenalty(2, 16);   % Gaussian constants + penalty, n = 16, k = beta + variance
                ok = isreal(A) && isfinite(A(1, 2)) && A(1, 5) == 2 && abs(A(1, 4) - part2) < 1e-6 && ...
                     L(1, 2) == dval && L(1, 3) == dval && L(1, 4) > 0;
                if dval == 1, ok = ok && L(1, 4) == 1; end
                detail = sprintf('AICc = %.4g, numparams = %d (expected 2), part2 = %.4f (expected %.4f); noiseLog = [%s]', ...
                    A(1, 2), A(1, 5), A(1, 4), part2, num2str(L(1, :), '%.4g '));

            case 6
                % Step 1: least-squares fit of the real data, to get a smooth model curve.
                qdfQuiet(@() qdfTestOptions(fit0, 'options_test_u6a', 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', 'model.name = ''SEIRu6a'';'));
                rng(1, 'twister');
                [~, ~, ~, curvesFit] = evalc('Run_Fit_ODEModel(@options_test_u6a, 1, 1, 17)');
                obs = load(fullfile('input', 'curve-flu1918SF.txt'));
                smooth = round(median(curvesFit(1:17, :), 2));
                smooth(1) = obs(1, 2);                            % keep the first observation
                smooth = max(smooth, 1);
                writematrix([obs(1:17, 1) smooth], fullfile('input', 'curve-test-smooth.txt'), 'Delimiter', 'tab');
                % Step 2: fit that smooth series with NB bootstrap noise requested.
                qdfQuiet(@() qdfTestOptions(fit0, 'options_test_u6', 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', 'model.name = ''SEIRu6'';', ...
                    dLine, 'dist1 = 2; % Error structure type', ...
                    'cadfilename1 = ''curve-flu1918SF'';', 'cadfilename1 = ''curve-test-smooth'';'));
                rng(1, 'twister');
                qdfQuiet(@() Run_Fit_ODEModel(@options_test_u6, 1, 1, 17));
                qdfQuiet(@() plotFit_ODEModel(@options_test_u6, 1, 1, 17));   % must find its files
                close all
                L = qdfNoiseLog('bootstraps-ODEModel-*SEIRu6-*.mat');
                nDist2 = numel(dir(fullfile('output', '*SEIRu6-*-method-0-dist-2-*')));
                nDist1 = numel(dir(fullfile('output', '*SEIRu6-*-method-0-dist-1-*')));
                if L(1, 6) > 1
                    ok = false;
                    detail = sprintf(['the fallback was not triggered: the dispersion of the smooth series is %.3g (> 1), ' ...
                        'so this run does not test it; noiseLog = [%s]'], L(1, 6), num2str(L(1, :), '%.4g '));
                else
                    ok = L(1, 2) == 2 && L(1, 3) == 1 && L(1, 4) == 1 && nDist2 > 0 && nDist1 == 0;
                    detail = sprintf('noiseLog = [%s]; %d files named -dist-2-, %d named -dist-1-; plotFit found its files', ...
                        num2str(L(1, :), '%.4g '), nDist2, nDist1);
                end

            case 7
                qdfQuiet(@() qdfTestOptions(fit3, 'options_test_u7', 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', 'model.name = ''SEIRu7'';', ...
                    'params.fixI0 = 1;', 'params.fixI0 = 0;'));
                rng(1, 'twister');
                [~, A] = evalc('Run_Fit_ODEModel(@options_test_u7, 1, 1, 17)');
                pen = qdfPenalty(3, 17);                          % beta, initial value, alpha; all 17 points count
                ok = isreal(A) && isfinite(A(1, 2)) && A(1, 5) == 3 && abs(A(1, 4) - pen) < 1e-6;
                detail = sprintf('AICc = %.4g, numparams = %d (expected 3), penalty = %.4f (expected %.4f)', ...
                    A(1, 2), A(1, 5), A(1, 4), pen);

            case 8
                obs = load(fullfile('input', 'curve-flu1918SF.txt'));
                writematrix([obs(:, 1) obs(:, 2) obs(:, 2)], fullfile('input', 'curve-test-two.txt'), 'Delimiter', 'tab');
                qdfQuiet(@() qdfTestOptions(fit3, 'options_test_u8', 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', 'model.name = ''SEIRu8'';', ...
                    'cadfilename1 = ''curve-flu1918SF'';', 'cadfilename1 = ''curve-test-two'';', ...
                    'vars.fit_index = 5;', 'vars.fit_index = [5 3];', ...
                    'vars.fit_diff = 1;', 'vars.fit_diff = [1 0];'));
                rng(1, 'twister');
                [~, A] = evalc('Run_Fit_ODEModel(@options_test_u8, 1, 1, 17)');
                pen = qdfPenalty(2, 32);                          % 2 x 17 points, the first of each series fixed
                n5 = numel(dir(fullfile('output', '*SEIRu8*-fit_index-5-*')));
                n3 = numel(dir(fullfile('output', '*SEIRu8*-fit_index-3-*')));
                ok = isreal(A) && isfinite(A(1, 2)) && A(1, 5) == 2 && abs(A(1, 4) - pen) < 1e-6 && n5 > 0 && n3 > 0;
                detail = sprintf('AICc = %.4g, numparams = %d (expected 2), penalty = %.4f (expected %.4f); files: %d for series 1, %d for series 2', ...
                    A(1, 2), A(1, 5), A(1, 4), pen, n5, n3);

            case 9
                qdfQuiet(@() qdfTestOptions(fc3, 'options_test_u9', 'B=300;', 'B=20;', ...
                    'model.name=''SEIR model'';', 'model.name=''SEIRu9'';'));
                rng(1, 'twister');
                [~, A] = evalc('Run_Forecasting_ODEModel(@options_test_u9, 1, 3, 17, 10)');
                L = qdfNoiseLog('bootstraps-ODEModel-*SEIRu9*.mat');
                d = dir(fullfile('output', 'parameters-ODEModel-*SEIRu9*.mat'));
                [~, k] = max([d.datenum]);
                S = load(fullfile(d(k).folder, d(k).name), 'PII_rangeFrac', 'PII_scale', 'param_estims');
                qdfQuiet(@() plotForecast_ODEModel(@options_test_u9, 1, 3, 17, 10));
                close all
                ok = size(A, 1) == 3 && isequal(A(:, 1), (1:3)') && all(isfinite(A(:, 2))) && ...
                     size(L, 1) == 3 && isequal(L(:, 1), (1:3)') && ...
                     size(S.PII_rangeFrac, 1) == 3 && size(S.PII_scale, 1) == 3 && size(S.param_estims, 3) == 3;
                detail = sprintf('AICc per window = [%s]; noiseLog rows = %d; PII rows = %d; plotForecast ran for all windows', ...
                    num2str(reshape(A(:, 2), 1, []), '%.4g '), size(L, 1), size(S.PII_rangeFrac, 1));

            case 10
                qdfQuiet(@() qdfTestOptions(pii3, 'options_test_u10', 'B=300;', 'B=20;'));
                if exist(outPII, 'dir') == 7, rmdir(outPII, 's'); end
                clear fit_model                                    % reset the once-per-session logscale warning
                txt = evalc('Run_PracticalIdentifiability_ODEModel(@options_test_u10, 30, 5, 2, [], outPII)');
                nWarn = numel(strfind(txt, 'searched on the LINEAR scale'));
                [~, S] = evalc('summarizeIdentifiability(outPII, 30, 2, [1 3], [0.18 5], {''r'', ''alpha''})');
                cols = S.Properties.VariableNames;
                need = {'median_PII', 'n_PII_undefined', 'median_rangeFrac', 'PII_scale', 'coverage95'};
                ok = istable(S) && numel(cols) == 17 && all(ismember(need, cols)) && height(S) == 2 && ...
                     all(S.files == 2) && S.PII_scale(1) == 0 && S.PII_scale(2) == 1 && ...
                     all(isfinite(S.median_rangeFrac)) && isfile(fullfile(outPII, 'recovery_summary.csv'));
                detail = sprintf(['%d columns, %d rows, files = [%s]; r: scale %g, range fraction %.4g, PII %.4g, undefined %d; ' ...
                    'alpha: scale %g, PII %.4g; logscale warning printed %d time(s)'], ...
                    numel(cols), height(S), num2str(reshape(S.files, 1, []), '%d '), S.PII_scale(1), S.median_rangeFrac(1), S.median_PII(1), ...
                    S.n_PII_undefined(1), S.PII_scale(2), S.median_PII(2), nWarn);
        end
    catch ME
        ok = false;
        where = '';
        if ~isempty(ME.stack)
            where = sprintf(' [%s, line %d]', ME.stack(1).name, ME.stack(1).line);
        end
        detail = ['ERROR: ' qdfFirstLine(ME.message) where];
        close all
    end
    close all
    if ok, verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('%s (%.0f s)\n', verdict, toc(t0));
    results(end+1, :) = {t, verdict, detail}; %#ok<AGROW>
end

fprintf('\n=============== QuantDiffForecast smoke test 2: summary ===============\n');
for r = 1:size(results, 1)
    fprintf('%s  %d  %s\n      %s\n', results{r, 2}, results{r, 1}, names{results{r, 1}}, results{r, 3});
end
fprintf('%d of %d tests passed.\n', sum(strcmp(results(:, 2), 'PASS')), size(results, 1));

if cleanup
    optFiles = {'options_test_u1.m', 'options_test_u2.m', 'options_test_u3.m', 'options_test_u4.m', ...
        'options_test_u5.m', 'options_test_u6.m', 'options_test_u6a.m', 'options_test_u7.m', ...
        'options_test_u8.m', 'options_test_u9.m', 'options_test_u10.m'};
    for k = 1:numel(optFiles)
        if exist(fullfile(pwd, optFiles{k}), 'file') == 2, delete(optFiles{k}); end
    end
    inFiles = {'curve-test-smooth.txt', 'curve-test-two.txt'};
    for k = 1:numel(inFiles)
        f = fullfile('input', inFiles{k});
        if exist(f, 'file') == 2, delete(f); end
    end
    d = dir(fullfile('output', '*SEIRu*'));
    for k = 1:numel(d)
        if ~d(k).isdir, delete(fullfile(d(k).folder, d(k).name)); end
    end
    if exist(outPII, 'dir') == 7, rmdir(outPII, 's'); end
    rehash
    fprintf('Test options files and their outputs were deleted.\n');
end
end


% ------------------------------ local helpers ------------------------------

function qdfQuiet(f)
% Run F with its console output hidden.
evalc('f();');
end

function p = qdfPenalty(k, n)
% AICc penalty for k estimated quantities and n observations.
p = 2 * k + (2 * k * (k + 1)) / (n - k - 1);
end

function L = qdfNoiseLog(pattern)
% noiseLog saved in the newest bootstraps file matching PATTERN.
d = dir(fullfile('output', pattern));
if isempty(d)
    error('qdfSmokeTest2:noFile', 'No file in output matches %s.', pattern);
end
[~, k] = max([d.datenum]);
S = load(fullfile(d(k).folder, d(k).name), 'noiseLog');
if ~isfield(S, 'noiseLog')
    error('qdfSmokeTest2:noNoiseLog', '%s has no noiseLog variable.', d(k).name);
end
L = S.noiseLog;
end

function s = qdfFirstLine(msg)
% First line of a (possibly multi-line) message, shortened.
parts = regexp(msg, '\n', 'split');
s = strtrim(parts{1});
if numel(s) > 200, s = [s(1:200) ' ...']; end
end
