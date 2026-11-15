function results = qdfSmokeTest3(tests, cleanup)
%QDFSMOKETEST3  Third set of MATLAB checks: the identifiability plots and two harder models.
%
%   qdfSmokeTest3                 run all four tests, then delete the test files
%   qdfSmokeTest3(1)              run only test 1
%   qdfSmokeTest3(1:4, false)     run all tests and keep the test files and outputs
%   results = qdfSmokeTest3(...)  also return the summary as a cell array
%
%   Run it from the toolbox code folder. It needs qdfTestOptions.m in the same
%   folder. Each test uses quick copies (20 bootstrap fits, own model names), so
%   your real results in 'output' are not touched.
%
%   Tests
%     1  plotPracticalIdentifiabilityResults: a copy of the script is run with
%        2 replicates of the exponential-growth study (window 30, serial); the
%        result files, recovery_summary.csv and the PDF figures must appear
%     2  weakly identified model: the SEIR example with beta, kappa AND gamma
%        estimated from 17 points. The run must complete and plot; the summary
%        line reports the identifiability diagnostics
%     3  profile likelihood of kappa in that model (expected to be flat or
%        bound-limited, not a finite interval)
%     4  a different dataset and a larger model: COVID-19 cases in Spain,
%        SEIR with under-reporting, five estimated parameters plus the
%        negative binomial dispersion, 60-point window
%
%   Tests 2 and 4 PASS when the run and the plot complete and the counts are
%   right. The diagnostics in their detail lines (unresolved bootstrap fits,
%   bound hits, PII ratio) describe how well the model is identified; they are
%   for reading, not part of the verdict.
%
%   The console output of the runners is hidden. A FAIL line carries the error
%   message or the value that was not as expected: send those lines back.

if nargin < 1 || isempty(tests),   tests = 1:4;     end
if nargin < 2 || isempty(cleanup), cleanup = true;  end

if exist('qdfTestOptions', 'file') ~= 2
    error('qdfSmokeTest3:helper', 'qdfTestOptions.m must be in the current folder or on the path.');
end
if exist('Run_Fit_ODEModel', 'file') ~= 2
    error('qdfSmokeTest3:folder', 'Run this from the toolbox code folder (Run_Fit_ODEModel.m not found).');
end
if ~isnumeric(tests) || any(tests(:) < 1) || any(tests(:) > 4) || any(tests(:) ~= fix(tests(:)))
    error('qdfSmokeTest3:tests', 'tests must be whole numbers between 1 and 4.');
end

names = { ...
    'plotPracticalIdentifiabilityResults (2 replicates)', ...
    'weakly identified SEIR (beta, kappa, gamma free)', ...
    'profile likelihood of kappa in that model', ...
    'COVID-19 Spain, 5 parameters + dispersion, 60 points'};
results = cell(0, 3);
fit3  = 'options_fit_SEIR_flu1918_dist1_3';
pii3  = 'options_forecast_PII_EXPO_r_dist1_3';
covid = 'options_fit_SEIR_unreported_covid_spain_dist1_3';
weakArgs = {'B = 300;', 'B = 20;', 'model.name = ''SEIR model'';', 'model.name = ''SEIRv2'';', ...
            'params.fixed = [0, 1, 1, 1];', 'params.fixed = [0, 0, 0, 1];'};

for t = tests(:)'
    fprintf('Test %d: %s ... ', t, names{t});
    t0 = tic;
    try
        switch t
            case 1
                % Options copy with its own model name, so the study gets its own folder.
                qdfQuiet(@() qdfTestOptions(pii3, 'options_test_v1', 'B=300;', 'B=20;', ...
                    'model.name=''EXP-PII-r-dist1-3'';', 'model.name=''EXPtest-PII-r-dist1-3'';'));
                % Copy of the script with small settings; the leading 'clear' is removed.
                txt = fileread('plotPracticalIdentifiabilityResults.m');
                txt = qdfReplace(txt, 'clear; close all; clc;', 'close all;');
                txt = qdfReplace(txt, 'options_handle = @options_forecast_PII_EXPO_r_dist1_3;', 'options_handle = @options_test_v1;');
                txt = qdfReplace(txt, 'windowsize1 = 30:10:50;', 'windowsize1 = 30;');
                txt = qdfReplace(txt, 'num_replicates = 20;', 'num_replicates = 2;');
                txt = qdfReplace(txt, 'use_parallel = true;', 'use_parallel = false;');
                fid = fopen('qdf_test_plotPII.m', 'w');  fwrite(fid, txt);  fclose(fid);
                rehash
                old = dir(fullfile('output', 'EXPtest_*'));
                for k = 1:numel(old)
                    if old(k).isdir, rmdir(fullfile(old(k).folder, old(k).name), 's'); end
                end
                qdfRunScript('qdf_test_plotPII');
                labels = qdfAxisLabels();
                close all
                d = dir(fullfile('output', 'EXPtest_*'));
                d = d([d.isdir]);
                if isempty(d)
                    error('qdfSmokeTest3:noFolder', 'The script did not create an output folder named EXPtest_*.');
                end
                folder = fullfile(d(1).folder, d(1).name);
                nRes = numel(dir(fullfile(folder, 'results_W30_rep*.mat')));
                pdfs = dir(fullfile(folder, '*.pdf'));
                hasGrid = any(strncmp({pdfs.name}, 'CI_grid_param_', 14));
                hasPII  = any(~strncmp({pdfs.name}, 'CI_grid_param_', 14));
                csv = fullfile(folder, 'recovery_summary.csv');
                nCols = 0;
                if isfile(csv)
                    T = readtable(csv);  nCols = width(T);
                end
                ok = nRes == 2 && hasGrid && hasPII && nCols == 17;
                detail = sprintf('folder %s: %d result files, %d PDF figures (%s), recovery_summary.csv with %d columns; y-axis labels: %s', ...
                    d(1).name, nRes, numel(pdfs), strjoin({pdfs.name}, ', '), nCols, labels);

            case 2
                qdfQuiet(@() qdfTestOptions(fit3, 'options_test_v2', weakArgs{:}));
                rng(1, 'twister');
                [txt, A] = evalc('Run_Fit_ODEModel(@options_test_v2, 1, 1, 17)');
                before = qdfFileVars('parameters-ODEModel-*SEIRv2*.mat');
                qdfQuiet(@() plotFit_ODEModel(@options_test_v2, 1, 1, 17));
                close all
                after = qdfFileVars('parameters-ODEModel-*SEIRv2*.mat');
                pen = qdfPenalty(4, 16);                           % beta, kappa, gamma, alpha; 16 informative points
                sci = qdfReadCsv('SCIs-rollingwindow-*SEIRv2*.csv');
                ok = isreal(A) && isfinite(A(1, 2)) && A(1, 5) == 4 && abs(A(1, 4) - pen) < 1e-6 && ...
                     isequal(before, after) && isreal(sci) && ~any(isinf(sci(:)));
                detail = sprintf('AICc = %.4g, numparams = %d (expected 4), penalty = %.4f (expected %.4f); %s', ...
                    A(1, 2), A(1, 5), A(1, 4), pen, qdfDiagnostics('*SEIRv2*', txt));

            case 3
                if exist('options_test_v2', 'file') ~= 2
                    qdfQuiet(@() qdfTestOptions(fit3, 'options_test_v2', weakArgs{:}));
                end
                rng(1, 'twister');
                [~, PL] = evalc(['Run_ProfileLikelihood_ODEModel(@options_test_v2, 1, 17, ' ...
                    'struct(''paramIndex'', 2, ''numstartpoints'', 3, ''printscreen'', 0))']);
                st = PL.param(1).status;
                ok = ischar(st) && ~isempty(st) && numel(PL.param(1).CI) == 2;
                detail = sprintf('reference: %s; status: %s; CI = [%.4g, %.4g]', ...
                    char(PL.referenceStatus), st, PL.param(1).CI(1), PL.param(1).CI(2));

            case 4
                qdfQuiet(@() qdfTestOptions(covid, 'options_test_v4', 'B=300;', 'B=20;', ...
                    'model.name=''SEIR covid model with underreporting'';', 'model.name=''COVIDv4'';'));
                clear fit_model                                    % reset the once-per-session logscale warning
                rng(1, 'twister');
                [txt, A] = evalc('Run_Fit_ODEModel(@options_test_v4, 1, 1, 60)');
                nWarn = numel(strfind(txt, 'searched on the LINEAR scale'));
                before = qdfFileVars('parameters-ODEModel-*COVIDv4*.mat');
                qdfQuiet(@() plotFit_ODEModel(@options_test_v4, 1, 1, 60));
                close all
                after = qdfFileVars('parameters-ODEModel-*COVIDv4*.mat');
                d = dir(fullfile('output', 'parameters-ODEModel-*COVIDv4*.mat'));
                [~, k] = max([d.datenum]);
                S = load(fullfile(d(k).folder, d(k).name), 'PII_scale');
                pen = qdfPenalty(6, 59);                           % 5 parameters + alpha; 60 points, the first fixed
                sci = qdfReadCsv('SCIs-rollingwindow-*COVIDv4*.csv');
                scaleOK = isequal(S.PII_scale(1, 1:5), [1 0 1 0 0]);   % lower bounds 0.001, 0, 0.5, 0, 0
                ok = isreal(A) && isfinite(A(1, 2)) && A(1, 5) == 6 && abs(A(1, 4) - pen) < 1e-6 && ...
                     isequal(before, after) && isreal(sci) && ~any(isinf(sci(:))) && scaleOK;
                detail = sprintf(['AICc = %.4g, numparams = %d (expected 6), penalty = %.4f (expected %.4f); ' ...
                    'PII scale of the 5 parameters = [%s] (expected 1 0 1 0 0); logscale warning printed %d time(s); %s'], ...
                    A(1, 2), A(1, 5), A(1, 4), pen, num2str(S.PII_scale(1, 1:5), '%d '), nWarn, qdfDiagnostics('*COVIDv4*', txt));
        end
    catch ME
        ok = false;
        where = '';
        if ~isempty(ME.stack)
            where = sprintf(' [%s, line %d]', ME.stack(1).name, ME.stack(1).line);
        end
        detail = ['ERROR: ' qdfFirstLine(ME.message) where];
    end
    close all
    if ok, verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('%s (%.0f s)\n', verdict, toc(t0));
    results(end+1, :) = {t, verdict, detail}; %#ok<AGROW>
end

fprintf('\n=============== QuantDiffForecast smoke test 3: summary ===============\n');
for r = 1:size(results, 1)
    fprintf('%s  %d  %s\n      %s\n', results{r, 2}, results{r, 1}, names{results{r, 1}}, results{r, 3});
end
fprintf('%d of %d tests passed.\n', sum(strcmp(results(:, 2), 'PASS')), size(results, 1));

if cleanup
    files = {'options_test_v1.m', 'options_test_v2.m', 'options_test_v4.m', 'qdf_test_plotPII.m'};
    for k = 1:numel(files)
        if exist(fullfile(pwd, files{k}), 'file') == 2, delete(files{k}); end
    end
    pats = {'*SEIRv2*', '*COVIDv4*'};
    for p = 1:numel(pats)
        d = dir(fullfile('output', pats{p}));
        for k = 1:numel(d)
            if ~d(k).isdir, delete(fullfile(d(k).folder, d(k).name)); end
        end
    end
    d = dir(fullfile('output', 'EXPtest_*'));
    for k = 1:numel(d)
        if d(k).isdir, rmdir(fullfile(d(k).folder, d(k).name), 's'); end
    end
    rehash
    fprintf('Test files and their outputs were deleted.\n');
end
end


% ------------------------------ local helpers ------------------------------

function qdfQuiet(f)
% Run F with its console output hidden.
evalc('f();');
end

function qdfRunScript(scriptName)
% Run a script in this function's own workspace (so a 'clear' or the script's
% variables cannot touch the caller), with its console output hidden.
evalc(scriptName);
end

function txt = qdfReplace(txt, oldText, newText)
% Replace OLDTEXT by NEWTEXT; error if OLDTEXT is not there.
if isempty(strfind(txt, oldText)) %#ok<STREMP>
    error('qdfSmokeTest3:notFound', 'Text not found in plotPracticalIdentifiabilityResults.m: "%s"', oldText);
end
txt = strrep(txt, oldText, newText);
end

function s = qdfAxisLabels()
% Y-axis labels of the open figures, as one line of text.
s = 'no open figure';
ax = findall(0, 'Type', 'axes');
if isempty(ax), return, end
labs = {};
for k = 1:numel(ax)
    try
        v = get(get(ax(k), 'YLabel'), 'String');
        if iscell(v), v = v{1}; end
        v = char(v);
        if ~isempty(v), labs{end+1} = v; end %#ok<AGROW>
    catch
        % ignore axes without a readable label
    end
end
if isempty(labs)
    s = 'none found';
else
    s = strjoin(unique(labs), ' | ');
end
end

function p = qdfPenalty(k, n)
% AICc penalty for k estimated quantities and n observations.
p = 2 * k + (2 * k * (k + 1)) / (n - k - 1);
end

function v = qdfFileVars(pattern)
% Sorted variable names of the newest file in 'output' matching PATTERN.
d = dir(fullfile('output', pattern));
if isempty(d)
    error('qdfSmokeTest3:noFile', 'No file in output matches %s.', pattern);
end
[~, k] = max([d.datenum]);
v = sort(who('-file', fullfile(d(k).folder, d(k).name)));
end

function x = qdfReadCsv(pattern)
% Numeric content of the newest CSV in 'output' matching PATTERN.
d = dir(fullfile('output', pattern));
if isempty(d)
    error('qdfSmokeTest3:noFile', 'No file in output matches %s.', pattern);
end
[~, k] = max([d.datenum]);
x = readmatrix(fullfile(d(k).folder, d(k).name));
end

function s = qdfDiagnostics(pattern, txt)
% One line describing how well the fit behind the files matching PATTERN is
% identified, from the runner's console text and its saved diagnostics.
unres = regexp(txt, 'Unresolved bootstrap datasets:\s*(\d+)', 'tokens', 'once');
reqd  = regexp(txt, 'Requested bootstrap datasets:\s*(\d+)', 'tokens', 'once');
if isempty(unres), unres = {'?'}; end
if isempty(reqd),  reqd  = {'?'}; end

d = dir(fullfile('output', ['parameters-ODEModel-' pattern '.mat']));
[~, k] = max([d.datenum]);
P = load(fullfile(d(k).folder, d(k).name));
hit = max([P.fracLB(:); P.fracUB(:)], [], 'omitnan');
ratio = max(P.PII_ratio(:), [], 'omitnan');
nOut = sum(P.estimOutsideCI(:) > 0);

conv = NaN;  agree = NaN;  grad = NaN;
if isfield(P, 'numericalAudit') && isfield(P.numericalAudit, 'originalFits') && ~isempty(P.numericalAudit.originalFits)
    F = P.numericalAudit.originalFits{1};
    if isfield(F, 'converged'),             conv  = double(F.converged);        end
    if isfield(F, 'startAgreement'),        agree = F.startAgreement;           end
    if isfield(F, 'scaledGradientInfNorm'), grad  = F.scaledGradientInfNorm;    end
end

s = sprintf(['diagnostics: unresolved bootstrap fits %s of %s; best fit converged = %g, start agreement = %.2f, ' ...
    'scaled gradient = %.2g; largest share of bootstrap fits on a bound = %.2f; largest PII ratio = %.2f; ' ...
    'estimates outside their own interval = %d'], unres{1}, reqd{1}, conv, agree, grad, hit, ratio, nOut);
end

function s = qdfFirstLine(msg)
% First line of a (possibly multi-line) message, shortened.
parts = regexp(msg, '\n', 'split');
s = strtrim(parts{1});
if numel(s) > 200, s = [s(1:200) ' ...']; end
end
