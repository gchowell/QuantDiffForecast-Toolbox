% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
%
function [P, residual, fitcurve, forecastcurve, timevect2,initialguess,fval,F1,F2,fitDiagnostics]=fit_model(data1,params0,numstartpoints,DT,modelX,paramsX,varsX,forecastingperiod)

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

% Validate input data before optimizing: non-finite values (e.g., in a
% bootstrap replicate curve) would make every local solver run fail
% silently inside MultiStart, so fail here with a clear message instead.
if any(~isfinite(data1(:)))
    error('fit_model:invalidData', ...
        ['fit_model received non-finite data (NaN/Inf) in data1. If this is a ', ...
        'bootstrap replicate, check the replicate-generation step ', ...
        '(AddErrorStructure) for invalid draws.']);
end

% Define the objective function handle
f = @parameterSearchODE;

% --- bound hardening: orientation, LB/UB consistency, seed projection ---
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

% Project seed z into [LB,UB] and repair NaNs
z0 = min(max(z,LB),UB);
mid = (LB + UB)/2;
nanZ = isnan(z0); if any(nanZ), z0(nanZ) = mid(nanZ); end

% Typical parameter magnitudes for problem scaling: the (projected) seed
% where it is nonzero, otherwise the bound midpoint, floored at 1e-3.
typx = max(abs(z0), 1e-3);
typx(z0 == 0) = max(abs(mid(z0 == 0)), 1e-3);

% Optimizer configuration: tight stopping tolerances (1e-6) with a large
% evaluation budget so runs from different start points converge to the
% same optimum; forward finite differences for gradients on an
% ODE-based objective; TypicalX/ScaleProblem so parameters of very
% different magnitudes (rates vs. population sizes) stay well
% conditioned.
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...                          % sqp, interior-point
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'FiniteDifferenceType', 'forward', ...
    'FiniteDifferenceStepSize',1e-4, ...
    'ScaleProblem', 'obj-and-constr', ...
    'TypicalX', typx, ...
    'MaxFunctionEvaluations', 20000, ...
    'MaxIterations', 20000);

% Define the optimization problem
problem = createOptimProblem('fmincon', 'objective', f, 'x0', z0, ...
                             'lb', LB, 'ub', UB, 'options', options);

%% === Start points: LHS over the box + jitters around the seed ===
nStarts = max(0, floor(numstartpoints));
span = UB - LB;

lhsOK = exist('lhsdesign','file') == 2;

% The full numstartpoints budget explores the box via LHS/maximin -
% box-wide coverage is what finds the global basin on multimodal
% likelihoods, so it is never reduced. Up to 3 extra starts are
% multiplicative jitters around the seed z0: for bootstrap refits the
% optimum lies near the seed by construction, so these warm starts
% converge quickly; for exploratory fits they add little cost.
nLHS = nStarts;
nJit = min(3, max(1, nStarts));

% One batch of start points, regenerated fresh on each retry attempt.
makeStarts = @() localMakeStarts(nLHS, nJit, d, LB, UB, span, z0, lhsOK);

starts = makeStarts();
initialguess = starts;
sp = CustomStartPointSet(starts);

% Setup MultiStart
ms = MultiStart('Display','off');
% Uncomment to run the independent local solves in parallel
% (requires the Parallel Computing Toolbox):


% If MultiStart reports failure (flagg < 0), retry with freshly drawn
% start points, up to maxAttempts passes. The best finite local minimum
% seen across attempts is tracked as a fallback.
maxAttempts = 3;
attempt = 0;
flagg = -Inf;
bestP = []; bestFval = Inf;

while flagg < 0 && attempt < maxAttempts

    attempt = attempt + 1;

    if attempt > 1
        starts = makeStarts();                 % fresh draws each retry
        initialguess = [initialguess; starts]; %#ok<AGROW>
        sp = CustomStartPointSet(starts);
    end

    [P, fval, flagg, outpt, allmins] = run(ms, problem, sp);

    % Track the best finite local minimum seen across attempts
    if ~isempty(allmins)
        [fv, idx] = min([allmins.Fval]);
        if isfinite(fv) && fv < bestFval
            bestFval = fv;
            bestP = allmins(idx).X;
        end
    end

end

if flagg < 0
    if ~isempty(bestP)
        warning('fit_model:multistart', ...
            'MultiStart did not converge cleanly after %d attempts; using best local minimum found (fval = %g).', ...
            maxAttempts, bestFval);
        P = bestP; fval = bestFval;
    else
        % On total failure, never return an empty P (downstream ODE solves would
        % crash inside the model right-hand side): fall back to the projected
        % seed z0 so all outputs have valid dimensions. MultiStart swallows
        % exceptions thrown by the objective, so probe the objective once
        % directly at the seed and report what happens.
        try
            ftest = f(z0);
            diagmsg = sprintf(['the objective evaluates to %g at the seed, so local runs ', ...
                'ended with nonpositive exit flags (no convergence).'], ftest);
        catch ME
            diagmsg = sprintf(['the objective THROWS an error at the seed: "%s". ', ...
                'Fix that underlying issue (often invalid data or model blow-up).'], ME.message);
        end
        warning('fit_model:multistart', ...
            ['MultiStart failed on all %d attempts and no local minimum was found; %s ', ...
            'Falling back to the seed parameters — results from this fit are NOT reliable.'], ...
            maxAttempts, diagmsg);
        P = z0;
        try, fval = f(z0); catch, fval = Inf; end
    end
end

% Finalize the returned parameters with the SAME evaluation as the objective.
% This replaces the previous separate tighter reconstruction; it does not add
% another optimization or another calibration solve. The returned fval is now
% the score of fitcurve under the existing observation/zero-handling rules.
optimizerFval = fval;
[fval, fitcurve, F1, IC] = parameterSearchODE(P);
[options_ode, numericalProfile] = quantdiffODEOptions(IC);
residual = fitcurve - ydata;

if nargout > 9
    fitDiagnostics = struct('numericalProfile',numericalProfile, ...
        'optimizerObjective',optimizerFval, 'finalObjective',fval, ...
        'objectiveDifference',fval-optimizerFval, ...
        'parameters',P, 'initialConditions',IC, ...
        'calibrationTimes',timevect, 'numberOfObservations',numel(ydata), ...
        'method',method1, 'zeroReplacement',0.001, ...
        'scoreAndCurveFromSameEvaluation',true, ...
        'optimizerExitflag',flagg, 'optimizationAttempts',attempt);
end

% Forecast handling
if forecastingperiod < 1
    forecastcurve = fitcurve;
    timevect2 = timevect;
    F2 = F1;
else
    % Build the forecast time grid in index units first, then scale by DT,
    % consistent with timevect = data1(:,1)*DT above; the first
    % length(timevect) points reproduce the calibration grid exactly and the
    % remainder extend it by forecastingperiod steps.
    timevect2 = (data1(1,1) : (data1(end,1) + forecastingperiod)) * DT;

    % Preserve the existing full-grid forecast reconstruction, now using the
    % common policy. Extending tspan can change the adaptive mesh: the
    % calibration portion of F2 need not be bitwise identical to F1. AICc is
    % calculated from F1/fitcurve, not from this extended-horizon solve.
    [tout2, F2] = ode15s(model.fc, timevect2, IC, options_ode, P, params.extra0);
    if size(F2,1) ~= numel(timevect2) || size(F2,2) ~= numel(IC) || ...
            numel(tout2) ~= numel(timevect2) || ...
            ~isequal(tout2(:),timevect2(:)) || ...
            ~isreal(F2) || any(~isfinite(F2(:)))
        error('QuantDiffForecast:InvalidForecastTrajectory', ...
            'The forecast solve must return complete finite real states at the requested times.');
    end
    forecastcurve = quantdiffObservationCurve(F2,vars);

end

end % main function


% <============================================================================>
% < Start-point batch generator:                                             ==>
% <  - nLHS points via LHS maximin over [LB,UB] (rand fallback)              ==>
% <  - nJit points as multiplicative log-space jitters around z0             ==>
% <  - always appends z0 itself, de-dups, clips to bounds, and guarantees    ==>
% <    at least one valid start                                              ==>
% <============================================================================>
function starts = localMakeStarts(nLHS, nJit, d, LB, UB, span, z0, lhsOK)

% Full-box exploration points
if nLHS > 0
    if lhsOK
        X = lhsdesign(nLHS, d, 'criterion','maximin','iterations',50);
        S1 = LB + X .* span;
    else
        S1 = LB + rand(nLHS, d) .* span;
    end
else
    S1 = zeros(0,d);
end

% Jittered warm starts around the seed: each coordinate is multiplied by
% 10^N(0,0.1) (i.e., roughly +/- 25% for one sigma), which respects the
% natural scale of each parameter. Coordinates at zero get a small additive
% perturbation instead.
if nJit > 0
    jitter = 10.^(0.1*randn(nJit, d));
    S2 = repmat(z0, nJit, 1) .* jitter;
    zeroCols = (z0 == 0);
    if any(zeroCols)
        S2(:, zeroCols) = repmat(z0(zeroCols), nJit, 1) + ...
            0.01*randn(nJit, sum(zeroCols)) .* repmat(span(zeroCols), nJit, 1);
    end
    S2 = min(max(S2, repmat(LB, nJit, 1)), repmat(UB, nJit, 1)); % clip to box
else
    S2 = zeros(0,d);
end

starts = [S1; S2; z0];

% Remove exact duplicate rows without changing parameter values.
starts = unique(starts, 'rows', 'stable');

% Keep rows inside bounds & finite
inB = all(starts >= (LB - 1e-12) & starts <= (UB + 1e-12), 2);
finiteReal = all(isfinite(starts), 2) & isreal(starts);
starts = starts(inB & finiteReal, :);

% Guarantee at least one valid start
if isempty(starts)
    starts = z0;
end

end
