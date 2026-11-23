function results = qdfTestInitialConditions()
%QDFTESTINITIALCONDITIONS  Checks of parameter-dependent initial conditions (vars.initialFcn).
%
%   qdfTestInitialConditions            run the ten checks and print a verdict each
%   results = qdfTestInitialConditions  also return {name, verdict, detail}
%
%   Run it from the toolbox code folder. It needs only base MATLAB (ode15s): no
%   optimization is run, so it takes a few seconds. The checks evaluate the
%   objective (parameterSearchODE) on rows 1-17 of the 1918 influenza series
%   with the SEIR example, Poisson likelihood, kappa = 1/1.9, gamma = 1/4.1.
%
%   Checks
%     1  no vars.initialFcn: initial conditions are vars.initial with the fitted
%        state set to X0 (the original behaviour)
%     2  a vars.initialFcn that returns vars.initial gives the same objective as
%        no function at all
%     3  reference value: the Poisson objective at beta = 0.76, N = 550000 agrees
%        with an independent Python/SciPy calculation (-836.7267)
%     4  without the function, (beta, N) and (2 beta, 2 N) give the same objective
%        (the redundancy the function removes)
%     5  with S(0) = N - I(0), the objective changes when N changes at fixed
%        beta/N, and the reference value at beta = 0.76/550, N = 1000 agrees with
%        Python (2670.0961)
%     6  the fitted state is set to X0 whatever the function returns for it
%     7  invalid results raise QuantDiffForecast:InvalidInitialFcn, and the
%        optimizer interface returns Inf for such a trial
%     8  vars.window is supplied by default and passed through; the rolling-window
%        warning fires only for unfitted states, a late window and no function
%     9  the example options_forecast_SEIR_flu1918_IC: its function reproduces
%        vars.initial at params.initial, keeps S + E + I + R = N, and starts a
%        window at row 21 with E(0) = 320/kappa = 608 and a population depleted
%        by the 837 earlier cases
%    10  objective values of that example agree with Python: the default function
%        at beta = 0.58, I_0 = 8 (-859.7854), independent of parameter 6; the
%        alternative with I_0 and E_0 as parameters at I_0 = 5.7, E_0 = 10.4
%        (-859.9134), and equal to the original SEIR example at I_0 = 4, E_0 = 0
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

global model params vars method1 timevect ydata %#ok<GVMIS>
saved = {model, params, vars, method1, timevect, ydata};
restore = onCleanup(@() localRestore(saved)); %#ok<NASGU>

basePath = fileparts(mfilename('fullpath'));
data = load(fullfile(basePath, 'input', 'curve-flu1918SF.txt'));
W = 17;

model = struct('fc', @SEIR1, 'name', 'SEIR model');
params = struct('label', {{'\beta', '\kappa', '\gamma', 'N'}}, 'num', 4, 'extra0', []);
vars0 = struct('label', {{'S', 'E', 'I', 'R', 'C'}}, 'initial', [550000 - 4, 0, 4, 0, 4], ...
    'fit_index', 5, 'fit_diff', 1);
method1 = 1;
timevect = data(1:W, 1);
ydata = data(1:W, 2);
y1 = ydata(1);

zOf = @(beta, N) [beta, 1/1.9, 1/4.1, N, y1, 1, 1];
coupled = @(p, x0, v) [p(4) - 4, 0, 4, 0, x0]; %#ok<INUSD>

results = cell(0, 3);

% 1 ---------------------------------------------------------------------------
vars = vars0;
IC = quantdiffInitialConditions(zOf(0.76, 550000), params, vars);
results = localAdd(results, '1 original behaviour', isequal(IC, [549996 0 4 0 y1]), mat2str(IC));

% 2 ---------------------------------------------------------------------------
vars = vars0;
fPlain = parameterSearchODE(zOf(0.76, 550000));
vars = vars0; vars.initialFcn = @(p, x0, v) v.initial;
fSame = parameterSearchODE(zOf(0.76, 550000));
results = localAdd(results, '2 identity function', fPlain == fSame, sprintf('%.10g vs %.10g', fPlain, fSame));

% 3 ---------------------------------------------------------------------------
ref = -836.7267;
results = localAdd(results, '3 reference value', abs(fPlain - ref) < 1e-3 * abs(ref), sprintf('%.6f (Python %.4f)', fPlain, ref));

% 4 ---------------------------------------------------------------------------
vars = vars0;
fScaled = parameterSearchODE(zOf(2 * 0.76, 2 * 550000));
results = localAdd(results, '4 redundancy without function', abs(fScaled - fPlain) < 1e-4 * abs(fPlain), ...
    sprintf('%.6f vs %.6f', fScaled, fPlain));

% 5 ---------------------------------------------------------------------------
vars = vars0; vars.initialFcn = coupled;
fBig = parameterSearchODE(zOf(0.76, 550000));
fSmall = parameterSearchODE(zOf(0.76 / 550, 1000));
ref5 = 2670.0961;
ok5 = abs(fBig - fPlain) < 1e-9 * abs(fPlain) && abs(fSmall - fBig) > 100 && abs(fSmall - ref5) < 1e-3 * abs(ref5);
results = localAdd(results, '5 coupling', ok5, sprintf('N=550000: %.6f, N=1000 at the same beta/N: %.6f (Python %.4f)', fBig, fSmall, ref5));

% 6 ---------------------------------------------------------------------------
vars = vars0; vars.initialFcn = @(p, x0, v) [p(4) - 4, 0, 4, 0, 999]; %#ok<INUSD>
IC = quantdiffInitialConditions(zOf(0.76, 550000), params, vars);
[~, curve] = parameterSearchODE(zOf(0.76, 550000));
results = localAdd(results, '6 fitted state pinned', IC(5) == y1 && curve(1) == y1, sprintf('IC(5) = %g, first fitted point = %g', IC(5), curve(1)));

% 7 ---------------------------------------------------------------------------
bad = {@(p, x0, v) [1 2 3], @(p, x0, v) [p(4) - 4, 0, -1, 0, x0], @(p, x0, v) [NaN, 0, 4, 0, x0], 'notAFunction'}; %#ok<INUSD>
nCaught = 0; nInf = 0;
for b = 1:numel(bad)
    vars = vars0; vars.initialFcn = bad{b};
    try
        quantdiffInitialConditions(zOf(0.76, 550000), params, vars);
    catch ME
        nCaught = nCaught + strcmp(ME.identifier, 'QuantDiffForecast:InvalidInitialFcn');
    end
    nInf = nInf + isinf(parameterSearchODE(zOf(0.76, 550000)));
end
vars = vars0; vars.initialFcn = coupled;
fInfeasible = parameterSearchODE(zOf(0.76, 2));          % N = 2 < I(0): S(0) would be negative
results = localAdd(results, '7 invalid results', nCaught == numel(bad) && nInf == numel(bad) && isinf(fInfeasible), ...
    sprintf('%d of %d errors, %d of %d Inf, infeasible N -> %g', nCaught, numel(bad), nInf, numel(bad), fInfeasible));

% 8 ---------------------------------------------------------------------------
vars = vars0; vars.initialFcn = @(p, x0, v) [p(4) - 4 - v.window.priorCum, 0, 4, v.window.priorCum, x0];
ICdefault = quantdiffInitialConditions(zOf(0.76, 550000), params, vars);
vars.window = struct('tstart', 21, 'rows', 21:37, 'time0', 20, 'firstObs', 320, 'priorCum', sum(data(1:20, 2)));
ICwin = quantdiffInitialConditions([0.76, 1/1.9, 1/4.1, 550000, 320, 1, 1], params, vars);
prior = sum(data(1:20, 2));
warnState = warning('off', 'QuantDiffForecast:rollingInitialStates');
f1 = quantdiffCheckInitialStates(vars0, 3);              % unfitted states, late window, no function
f2 = quantdiffCheckInitialStates(vars0, 1);              % first window only
f3 = quantdiffCheckInitialStates(vars, 3);               % function supplied
f4 = quantdiffCheckInitialStates(struct('initial', 5, 'fit_index', 1), 3);   % one-state model
warning(warnState);
ok8 = isequal(ICdefault, [549996 0 4 0 y1]) && isequal(ICwin, [550000 - 4 - prior, 0, 4, prior, 320]) && ...
    f1 && ~f2 && ~f3 && ~f4;
results = localAdd(results, '8 window context and warning', ok8, sprintf('prior cases = %g, S(0) = %g, flags = %d%d%d%d', prior, ICwin(1), f1, f2, f3, f4));

% 9 ---------------------------------------------------------------------------
% (the options function sets the global method1; it is set back to 1 for check 10)
[~, ~, ~, ~, ~, ~, modelE, paramsE, varsE] = options_forecast_SEIR_flu1918_IC();
paramsE.num = numel(paramsE.label);
Ntrue = paramsE.initial(4); kappaE = paramsE.initial(2); I0E = paramsE.initial(5);
ICe = quantdiffInitialConditions([paramsE.initial, varsE.initial(varsE.fit_index), 1, 1], paramsE, varsE);
varsW = varsE;
varsW.window = struct('tstart', 21, 'rows', 21:37, 'time0', 20, 'firstObs', 320, 'priorCum', prior);
ICe21 = quantdiffInitialConditions([paramsE.initial, 320, 1, 1], paramsE, varsW);
expect21 = [Ntrue - I0E - 320 / kappaE - prior, 320 / kappaE, I0E, prior, 320];
ok9 = max(abs(ICe - varsE.initial)) < 1e-9 && abs(sum(ICe(1:4)) - Ntrue) < 1e-9 && ...
    max(abs(ICe21 - expect21)) < 1e-9 && abs(sum(ICe21(1:4)) - Ntrue) < 1e-9 && ...
    ~quantdiffCheckInitialStates(varsE, 21);
results = localAdd(results, '9 example options file', ok9, sprintf('row 1: %s; row 21: %s', mat2str(ICe, 8), mat2str(ICe21, 8)));

% 10 --------------------------------------------------------------------------
model = modelE; params = paramsE; vars = varsE; method1 = 1;
fTied  = parameterSearchODE([0.58, 1/1.9, 1/4.1, 550000, 8, 4, y1, 1, 1]);
fTied2 = parameterSearchODE([0.58, 1/1.9, 1/4.1, 550000, 8, 50, y1, 1, 1]);     % parameter 6 is not read by the default function
refTied = -859.7854;
% the alternative kept as a comment in the options file: I_0 and E_0 both parameters
vars.initialFcn = @(p, x0, v) [p(4)-p(5)-p(6)-v.window.priorCum, p(6), p(5), v.window.priorCum, x0];
fNested = parameterSearchODE([0.76, 1/1.9, 1/4.1, 550000, 4, 0, y1, 1, 1]);     % I_0 = 4, E_0 = 0: the original SEIR example
fBoth   = parameterSearchODE([0.58, 1/1.9, 1/4.1, 550000, 5.7, 10.4, y1, 1, 1]);
refBoth = -859.9134;
ok10 = abs(fTied - refTied) < 1e-3 * abs(refTied) && fTied == fTied2 && ...
    abs(fNested - fPlain) < 1e-9 * abs(fPlain) && abs(fBoth - refBoth) < 1e-3 * abs(refBoth);
results = localAdd(results, '10 example objective values', ok10, ...
    sprintf('default: %.6f (Python %.4f); alternative: %.6f (Python %.4f); alternative at I_0=4, E_0=0: %.6f (original %.6f)', ...
    fTied, refTied, fBoth, refBoth, fNested, fPlain));

fprintf('\n%d of %d checks passed.\n', sum(strcmp(results(:, 2), 'PASS')), size(results, 1));
end


function results = localAdd(results, name, ok, detail)
if ok, verdict = 'PASS'; else, verdict = 'FAIL'; end
fprintf('%-34s %s   %s\n', name, verdict, detail);
results(end+1, :) = {name, verdict, detail};
end

function localRestore(saved)
global model params vars method1 timevect ydata %#ok<GVMIS>
[model, params, vars, method1, timevect, ydata] = saved{:};
end
