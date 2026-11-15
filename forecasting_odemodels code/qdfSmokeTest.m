function results = qdfSmokeTest(tests, cleanup)
%QDFSMOKETEST  Run the QuantDiffForecast checks that have not yet been exercised in MATLAB.
%
%   qdfSmokeTest                 run all eight tests, then delete the test files
%   qdfSmokeTest([1 3])          run only tests 1 and 3
%   qdfSmokeTest(1:8, false)     run all tests and keep the test files and outputs
%   results = qdfSmokeTest(...)  also return the summary as a cell array
%
%   Run it from the toolbox code folder. It needs qdfTestOptions.m in the same
%   folder. Each test uses a quick copy of a shipped options file (20 bootstrap
%   fits, its own model name), so your real results in 'output' are not touched.
%
%   Tests
%     1  plotFit_ODEModel keeps the variables the runner saved (-append save)
%     2  plotForecast_ODEModel keeps the variables the runner saved
%     3  least squares with negative binomial noise (method 0, dist1 2):
%        dispersion from residuals, file names keep '-dist-2-', plotFit finds them
%     4  SAD (method 6): Laplace scale from the mean absolute residual
%     5  profile likelihood, negative binomial fit
%     6  profile likelihood, least-squares fit (uses the effective sample size)
%     7  method1 = 2 stops at once with a clear message
%     8  negative counts with a count likelihood stop at once, listing the row
%
%   The console output of the runners is hidden; only a one-line verdict per
%   test and a summary are printed. A FAIL line carries the error message or the
%   value that was not as expected: send those lines back for diagnosis.

if nargin < 1 || isempty(tests),   tests = 1:8;     end
if nargin < 2 || isempty(cleanup), cleanup = true;  end

if exist('qdfTestOptions', 'file') ~= 2
    error('qdfSmokeTest:helper', 'qdfTestOptions.m must be in the current folder or on the path.');
end
if exist('Run_Fit_ODEModel', 'file') ~= 2
    error('qdfSmokeTest:folder', 'Run this from the toolbox code folder (Run_Fit_ODEModel.m not found).');
end
if ~isnumeric(tests) || any(tests(:) < 1) || any(tests(:) > 8) || any(tests(:) ~= fix(tests(:)))
    error('qdfSmokeTest:tests', 'tests must be whole numbers between 1 and 8.');
end

names = { ...
    'plotFit keeps the runner''s saved variables', ...
    'plotForecast keeps the runner''s saved variables', ...
    'least squares + NB noise (method 0, dist1 2)', ...
    'SAD (method 6)', ...
    'profile likelihood, negative binomial', ...
    'profile likelihood, least squares', ...
    'method1 = 2 stops early', ...
    'negative counts stop early'};
results = cell(0, 3);
fitSrc  = 'options_fit_SEIR_flu1918_dist1_3';
fitSrc0 = 'options_fit_SEIR_flu1918_dist1_0';
fcSrc   = 'options_forecast_SEIR_flu1918_dist1_3';

for t = tests(:)'
    fprintf('Test %d: %s ... ', t, names{t});
    t0 = tic;
    try
        switch t
            case 1
                qdfQuiet(@() qdfTestOptions(fitSrc, 'options_test_t1', 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', 'model.name = ''SEIRt1'';'));
                rng(1, 'twister');
                qdfQuiet(@() Run_Fit_ODEModel(@options_test_t1, 1, 1, 17));
                before = qdfFileVars('parameters-ODEModel-*SEIRt1*.mat');
                qdfQuiet(@() plotFit_ODEModel(@options_test_t1, 1, 1, 17));
                close all
                after = qdfFileVars('parameters-ODEModel-*SEIRt1*.mat');
                need = {'PII_rangeFrac', 'PII_scale', 'PII_ratio', 'numericalAudit', 'param_estims'};
                ok = isequal(before, after) && all(ismember(need, after));
                detail = sprintf('%d variables before plotting, %d after', numel(before), numel(after));
                sci = qdfReadCsv('SCIs-rollingwindow-*SEIRt1*.csv');
                if any(isinf(sci(:))) || ~isreal(sci)
                    ok = false; detail = [detail '; the SCIs CSV holds Inf or complex values'];
                end

            case 2
                qdfQuiet(@() qdfTestOptions(fcSrc, 'options_test_t2', 'B=300;', 'B=20;', ...
                    'model.name=''SEIR model'';', 'model.name=''SEIRt2'';'));
                rng(1, 'twister');
                qdfQuiet(@() Run_Forecasting_ODEModel(@options_test_t2, 1, 1, 17, 10));
                before = qdfFileVars('parameters-ODEModel-*SEIRt2*.mat');
                qdfQuiet(@() plotForecast_ODEModel(@options_test_t2, 1, 1, 17, 10));
                close all
                after = qdfFileVars('parameters-ODEModel-*SEIRt2*.mat');
                need = {'PII_rangeFrac', 'PII_scale', 'PII_ratio', 'numericalAudit', 'param_estims'};
                ok = isequal(before, after) && all(ismember(need, after));
                detail = sprintf('%d variables before plotting, %d after', numel(before), numel(after));

            case 3
                qdfQuiet(@() qdfTestOptions(fitSrc0, 'options_test_t3', 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', 'model.name = ''SEIRt3'';', ...
                    'dist1 = 0; % Error structure type', 'dist1 = 2; % Error structure type'));
                rng(1, 'twister');
                qdfQuiet(@() Run_Fit_ODEModel(@options_test_t3, 1, 1, 17));
                qdfQuiet(@() plotFit_ODEModel(@options_test_t3, 1, 1, 17));   % must find its files
                close all
                L = qdfNoiseLog('bootstraps-ODEModel-*SEIRt3*.mat');
                raw = L(1, 6);
                if isfinite(raw) && raw > 1
                    okNoise = L(1, 3) == 2 && abs(L(1, 4) - raw) < 1e-12;
                else
                    okNoise = L(1, 3) == 1 && L(1, 4) == 1;       % Poisson fallback
                end
                nDist2 = numel(dir(fullfile('output', '*SEIRt3*-method-0-dist-2-*')));
                nDist1 = numel(dir(fullfile('output', '*SEIRt3*-method-0-dist-1-*')));
                ok = L(1, 2) == 2 && okNoise && nDist2 > 0 && nDist1 == 0;
                detail = sprintf('noiseLog = [%s]; %d files named -dist-2-, %d named -dist-1-', ...
                    num2str(L(1, :), '%.4g '), nDist2, nDist1);

            case 4
                qdfQuiet(@() qdfTestOptions(fitSrc0, 'options_test_t4', 'B = 300;', 'B = 20;', ...
                    'model.name = ''SEIR model'';', 'model.name = ''SEIRt4'';', ...
                    'method1 = 0; % Parameter estimation method', 'method1 = 6; % Parameter estimation method'));
                rng(1, 'twister');
                [~, AICcs] = evalc('Run_Fit_ODEModel(@options_test_t4, 1, 1, 17)');
                L = qdfNoiseLog('bootstraps-ODEModel-*SEIRt4*.mat');
                ok = L(1, 2) == 6 && L(1, 3) == 6 && L(1, 4) > 0 && abs(L(1, 4) - L(1, 6)) < 1e-12 && ...
                     isfinite(AICcs(1, 2)) && isreal(AICcs(1, 2));
                detail = sprintf('noiseLog = [%s]; AICc = %.4g', num2str(L(1, :), '%.4g '), AICcs(1, 2));

            case 5
                if exist('options_test_t1', 'file') ~= 2
                    qdfQuiet(@() qdfTestOptions(fitSrc, 'options_test_t1', 'B = 300;', 'B = 20;', ...
                        'model.name = ''SEIR model'';', 'model.name = ''SEIRt1'';'));
                end
                rng(1, 'twister');
                [~, PL] = evalc(['Run_ProfileLikelihood_ODEModel(@options_test_t1, 1, 17, ' ...
                    'struct(''paramIndex'', 1, ''numstartpoints'', 3, ''printscreen'', 0))']);
                ok = ischar(PL.param(1).status) && ~isempty(PL.param(1).status) && numel(PL.param(1).CI) == 2;
                detail = sprintf('reference: %s; status: %s; CI = [%.4g, %.4g]', ...
                    char(PL.referenceStatus), PL.param(1).status, PL.param(1).CI(1), PL.param(1).CI(2));

            case 6
                if exist('options_test_t3', 'file') ~= 2
                    qdfQuiet(@() qdfTestOptions(fitSrc0, 'options_test_t3', 'B = 300;', 'B = 20;', ...
                        'model.name = ''SEIR model'';', 'model.name = ''SEIRt3'';', ...
                        'dist1 = 0; % Error structure type', 'dist1 = 2; % Error structure type'));
                end
                rng(1, 'twister');
                [~, PL] = evalc(['Run_ProfileLikelihood_ODEModel(@options_test_t3, 1, 17, ' ...
                    'struct(''paramIndex'', 1, ''numstartpoints'', 3, ''printscreen'', 0))']);
                ok = ischar(PL.param(1).status) && ~isempty(PL.param(1).status) && numel(PL.param(1).CI) == 2;
                detail = sprintf('reference: %s; status: %s; CI = [%.4g, %.4g]', ...
                    char(PL.referenceStatus), PL.param(1).status, PL.param(1).CI(1), PL.param(1).CI(2));

            case 7
                qdfQuiet(@() qdfTestOptions(fitSrc, 'options_test_t7', ...
                    'method1 = 3; % Parameter estimation method', 'method1 = 2; % Parameter estimation method'));
                [id, msg, secs] = qdfExpectError(@() Run_Fit_ODEModel(@options_test_t7, 1, 1, 17));
                ok = strcmp(id, 'QuantDiffForecast:UnsupportedMethod') && ~isempty(strfind(msg, 'method1 = 2 is not supported')); %#ok<STREMP>
                detail = sprintf('%s after %.1f s: %s', id, secs, qdfFirstLine(msg));

            case 8
                obs = load(fullfile('input', 'curve-flu1918SF.txt'));
                obs(5, 2) = -3;
                writematrix(obs, fullfile('input', 'curve-test-negative.txt'), 'Delimiter', 'tab');
                qdfQuiet(@() qdfTestOptions(fitSrc, 'options_test_t8', ...
                    'cadfilename1 = ''curve-flu1918SF'';', 'cadfilename1 = ''curve-test-negative'';'));
                [id, msg, secs] = qdfExpectError(@() Run_Fit_ODEModel(@options_test_t8, 1, 1, 17));
                ok = strcmp(id, 'QuantDiffForecast:NegativeCounts') && ~isempty(strfind(msg, 'row 5')) && ...
                     ~isempty(strfind(msg, '-3')); %#ok<STREMP>
                detail = sprintf('%s after %.1f s: %s', id, secs, qdfFirstLine(msg));

            otherwise
                error('qdfSmokeTest:test', 'There is no test %d (use 1 to 8).', t);
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
    if ok, verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('%s (%.0f s)\n', verdict, toc(t0));
    results(end+1, :) = {t, verdict, detail}; %#ok<AGROW>
end

fprintf('\n================ QuantDiffForecast smoke test: summary ================\n');
for r = 1:size(results, 1)
    fprintf('%s  %d  %s\n      %s\n', results{r, 2}, results{r, 1}, names{results{r, 1}}, results{r, 3});
end
fprintf('%d of %d tests passed.\n', sum(strcmp(results(:, 2), 'PASS')), size(results, 1));

if cleanup
    for t = [1 2 3 4 7 8]
        f = sprintf('options_test_t%d.m', t);
        if exist(fullfile(pwd, f), 'file') == 2, delete(f); end
    end
    f = fullfile('input', 'curve-test-negative.txt');
    if exist(f, 'file') == 2, delete(f); end
    for pat = {'*SEIRt1*', '*SEIRt2*', '*SEIRt3*', '*SEIRt4*'}
        d = dir(fullfile('output', pat{1}));
        for k = 1:numel(d)
            if ~d(k).isdir, delete(fullfile(d(k).folder, d(k).name)); end
        end
    end
    rehash
    fprintf('Test options files and their outputs were deleted.\n');
end
end


% ------------------------------ local helpers ------------------------------

function qdfQuiet(f)
% Run F with its console output hidden.
evalc('f();');
end

function v = qdfFileVars(pattern)
% Sorted variable names of the newest file in 'output' matching PATTERN.
d = dir(fullfile('output', pattern));
if isempty(d)
    error('qdfSmokeTest:noFile', 'No file in output matches %s.', pattern);
end
[~, k] = max([d.datenum]);
v = sort(who('-file', fullfile(d(k).folder, d(k).name)));
end

function L = qdfNoiseLog(pattern)
% noiseLog saved in the newest bootstraps file matching PATTERN.
d = dir(fullfile('output', pattern));
if isempty(d)
    error('qdfSmokeTest:noFile', 'No file in output matches %s.', pattern);
end
[~, k] = max([d.datenum]);
S = load(fullfile(d(k).folder, d(k).name), 'noiseLog');
if ~isfield(S, 'noiseLog')
    error('qdfSmokeTest:noNoiseLog', '%s has no noiseLog variable.', d(k).name);
end
L = S.noiseLog;
end

function x = qdfReadCsv(pattern)
% Numeric content of the newest CSV in 'output' matching PATTERN.
d = dir(fullfile('output', pattern));
if isempty(d)
    error('qdfSmokeTest:noFile', 'No file in output matches %s.', pattern);
end
[~, k] = max([d.datenum]);
x = readmatrix(fullfile(d(k).folder, d(k).name));
end

function [id, msg, secs] = qdfExpectError(f)
% Run F, which is expected to stop with an error; return its identifier,
% message and the time it took. No error -> id 'none'.
id = 'none'; msg = 'the call finished without an error';
t0 = tic;
try
    evalc('f();');
catch ME
    id = ME.identifier; msg = ME.message;
end
secs = toc(t0);
end

function s = qdfFirstLine(msg)
% First line of a (possibly multi-line) message, shortened.
parts = regexp(msg, '\n', 'split');
s = strtrim(parts{1});
if numel(s) > 200, s = [s(1:200) ' ...']; end
end
