% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [P, residual, fitcurve, forecastcurve, timevect2,initialguess,fval,F1,F2]=fit_model(data1,params0,numstartpoints,DT,modelX,paramsX,varsX,forecastingperiod)

global model params vars method1 timevect ydata

model=modelX;
params=paramsX;
vars=varsX;

% Calculate time vector from data column and a scaling factor DT
timevect = data1(:,1) * DT;

% Extract initial conditions from the second column onwards of the first row
I0 = data1(1, 2:end);

% Copy initial parameter values
z = params0;

% Adjust lower and upper bounds for parameters that are fixed
for i = 1:params.num
    if params.fixed(i)
        params.LB(i) = params0(i);
        params.UB(i) = params0(i);
    end
end

% Set extended bounds based on the method specified
switch method1
    case {0, 1,6}
        % Cases 0 and 1 have no extension on bounds
        LBe = [0 0];
        UBe = [0 0];
    case {3, 4}
        % Cases 3 and 4 allow wide variation
        LBe = [1e-8, 1];
        UBe = [1e4, 1];
    case 5
        % Case 5 has specific limits for d, assuming it should be >=1
        LBe = [1e-8, 0.6];
        UBe = [1e4, 1e3];
end

% Configure bounds based on whether initial conditions are fixed
if params.fixI0 == 1
    LB = [params.LB, I0, LBe];
    UB = [params.UB, I0, UBe];
else
    % If not fixed, use zero and sum of absolute values for initial conditions
    LB = [params.LB, zeros(1, length(I0))+0.001, LBe];
    UB = [params.UB, sum(abs(data1(:, 2:end))), UBe];
end

% ydata for fitting (vectorize if multiple series)
ydata = data1(:,2:end);
if length(I0)>1
    ydata = ydata(:);
end

% Define the objective function handle
f = @parameterSearchODE;

% Set optimization options for fmincon
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...                          % sqp, interior-point
    'StepTolerance', 1e-4, ...
    'FunctionTolerance', 1e-4, ...
    'OptimalityTolerance', 1e-4, ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 10000);

% Define the optimization problem
problem = createOptimProblem('fmincon', 'objective', f, 'x0', z, ...
                             'lb', LB, 'ub', UB, 'options', options);

%% === Smarter start points (LHS maximin + de-dup + include z; hardened) ===
d = numel(z);
LB = LB(:)'; UB = UB(:)'; z = z(:)';
if numel(LB) ~= d || numel(UB) ~= d
    error('Length mismatch: numel(z)=%d, numel(LB)=%d, numel(UB)=%d', d, numel(LB), numel(UB));
end
flipMask = LB > UB;
if any(flipMask)
    tmp = LB(flipMask); LB(flipMask) = UB(flipMask); UB(flipMask) = tmp;
end
if ~all(isfinite(LB)) || ~all(isfinite(UB))
    error('LB/UB must be finite.');
end

nStarts = max(0, floor(numstartpoints));
span = UB - LB;

% Project seed z into [LB,UB] and repair NaNs
z0 = min(max(z,LB),UB);
mid = (LB + UB)/2;
nanZ = isnan(z0); if any(nanZ), z0(nanZ) = mid(nanZ); end

lhsOK = exist('lhsdesign','file') == 2;
if lhsOK && nStarts > 0
    X = lhsdesign(nStarts, d, 'criterion','maximin','iterations',50);
    starts = LB + X .* span;
elseif nStarts > 0
    starts = LB + rand(nStarts, d) .* span;
else
    starts = zeros(0,d);
end

starts = [starts; z0];
starts = unique(round(starts,6), 'rows');

% Keep rows inside bounds & finite
inB = all(starts >= (LB - 1e-12) & starts <= (UB + 1e-12), 2);
finiteReal = all(isfinite(starts), 2) & isreal(starts);
starts = starts(inB & finiteReal, :);

% Guarantee at least one valid start
if isempty(starts)
    starts = z0;
end

initialguess = starts;
sp = CustomStartPointSet(starts);

% Setup MultiStart (quiet; parallel if pool exists)
useParallel = ~isempty(gcp('nocreate'));
useParallel = 0;
ms = MultiStart('Display','off', ...
                'UseParallel',useParallel, ...
                'StartPointsToRun','bounds-ineqs');

% Run MultiStart once with smarter starts
[P, fval, flagg, outpt, allmins] = run(ms, problem, sp);

% Initialize the options for the ODE solver (if any specific options needed, define here)
options_ode = [];

% Set initial conditions based on params.fixI0 flag
IC = vars.initial;
if params.fixI0 == 1
    IC(vars.fit_index) = I0;  % Fix the initial conditions to I0 for specified indices
else
    % If not fixed, use parameter values following the first 'num' parameters
    IC(vars.fit_index) = P(params.num + 1 : params.num + length(I0));
end

% Solve the differential equations using ode15s
[~, F] = ode15s(model.fc, timevect, IC, options_ode, P, params.extra0);
F1 = F;

% Build fitted curve (handles levels vs. diffs per variable)
yfit = zeros(length(ydata), 1);
currentEnd = 0;
for j = 1:length(vars.fit_index)
    if vars.fit_diff(j) == 1
        fitcurve = abs([F(1, vars.fit_index(j)); diff(F(:, vars.fit_index(j)))]);
    else
        fitcurve = F(:, vars.fit_index(j));
    end
    yfit(currentEnd + 1 : currentEnd + length(fitcurve)) = fitcurve;
    currentEnd = currentEnd + length(fitcurve);
end
fitcurve = yfit;

% Residuals
residual = fitcurve - ydata;

% Forecast handling
if forecastingperiod < 1
    forecastcurve = residual + ydata;
    timevect2 = timevect;
    F2 = F1;
else
    timevect2 = data1(1,1) : DT : (data1(end,1) + forecastingperiod * DT);
    [~, F2] = ode15s(model.fc, timevect2, IC, options_ode, P, params.extra0);

    yforecast = zeros(length(vars.fit_index) * length(timevect2), 1);
    currentEnd = 0;
    for j = 1:length(vars.fit_index)
        if vars.fit_diff(j) == 1
            forecastcurve = abs([F2(1, vars.fit_index(j)); diff(F2(:, vars.fit_index(j)))]);
        else
            forecastcurve = F2(:, vars.fit_index(j));
        end
        yforecast(currentEnd + 1 : currentEnd + length(forecastcurve)) = forecastcurve;
        currentEnd = currentEnd + length(forecastcurve);
    end
    forecastcurve = yforecast;
end
