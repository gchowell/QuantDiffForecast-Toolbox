% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
%
function [P, residual, fitcurve, forecastcurve, timevect2,initialguess,fval,F1,F2,fitDiagnostics]=fit_model(data1,params0,numstartpoints,DT,modelX,paramsX,varsX,forecastingperiod,pin)
% Optional 9th input PIN (added for Run_ProfileLikelihood_ODEModel): a vector the
% length of the full optimizer vector [params, X0, alpha, d]. Non-NaN entries pin
% that coordinate at the given value (LB = UB = value) for this call only, which
% lets any coordinate (including X0 and alpha) be held fixed while the rest are
% re-optimized. Omit or pass [] for the original behaviour; nothing else changes.
%
% Optional fields of paramsX (both default to the original behaviour):
%   paramsX.logscale       true  -> the optimizer works in log10 of every coordinate
%                                   whose lower bound is positive (rates, N, X0, alpha),
%                                   which conditions problems whose parameters span
%                                   several orders of magnitude. Bounds, start points
%                                   and all outputs stay in the original scale.
%   paramsX.scaleObjective false -> fmincon's 'ScaleProblem' is 'none' instead of
%                                   'obj-and-constr'. With 'obj-and-constr' the
%                                   objective is normalized by its starting value, so
%                                   for large |log-likelihoods| (long series, large
%                                   counts) a change of one deviance unit falls below
%                                   the 1e-6 tolerances and SQP stops at the seed.
%                                   'none' keeps the tolerances in objective units.
%   paramsX.polish         true  -> after MultiStart, run a derivative-free Nelder-Mead
%                                   polish (fminsearch) on the free coordinates from the
%                                   SQP solution and keep it if it lowers the objective.
%                                   Robust when |log L| is large and SQP's finite-
%                                   difference steps stall; costs a few hundred extra
%                                   objective evaluations per fit.

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
        UBe = [1e4, 5];
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

% Apply the optional PIN vector (profile-likelihood support): pinned coordinates
% get LB = UB = pin value. Applied after the fixed-parameter and X0 rules above,
% so it also overrides params.fixI0 for the pinned X0 slot.
if nargin >= 9 && ~isempty(pin)
    pin = pin(:)';
    if numel(pin) ~= numel(LB)
        error('fit_model:pinLength', 'pin must have %d entries (params, X0, alpha, d); got %d.', numel(LB), numel(pin));
    end
    pinned = ~isnan(pin);
    LB(pinned) = pin(pinned);
    UB(pinned) = pin(pinned);
    z(pinned)  = pin(pinned);
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
useLog = isfield(params, 'logscale') && ~isempty(params.logscale) && params.logscale;
if isfield(params, 'scaleObjective') && ~isempty(params.scaleObjective) && ~params.scaleObjective
    scaleProblem = 'none';
else
    scaleProblem = 'obj-and-constr';
end

% Optional log10 reparametrization of the optimizer's coordinates (see header).
% tf marks the transformed coordinates; LBq/UBq/q0 are the optimizer's box and seed.
if useLog
    tf = (LB > 0) & (UB > LB);
else
    tf = false(1, d);
end
toQ   = @(zz) localToQ(zz, tf);
fromQ = @(qq) localFromQ(qq, tf);
LBq = toQ(LB);  UBq = toQ(UB);  q0 = toQ(z0);
fq  = @(qq) f(fromQ(qq));
typq = max(abs(q0), 1e-3);
typq(q0 == 0) = max(abs((LBq(q0 == 0) + UBq(q0 == 0)) / 2), 1e-3);
if ~useLog, typq = typx; end

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...                          % sqp, interior-point
    'StepTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'FiniteDifferenceType', 'forward', ...
    'FiniteDifferenceStepSize',1e-4, ...
    'ScaleProblem', scaleProblem, ...
    'TypicalX', typq, ...
    'MaxFunctionEvaluations', 20000, ...
    'MaxIterations', 20000);

% Define the optimization problem (in the optimizer's coordinates)
problem = createOptimProblem('fmincon', 'objective', fq, 'x0', q0, ...
                             'lb', LBq, 'ub', UBq, 'options', options);

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

% One batch of start points, regenerated fresh on each retry attempt. With
% logscale on, the LHS points and jitters are drawn in the optimizer's log10
% coordinates (uniform over decades, additive jitter) and mapped back.
if useLog
    makeStarts = @() fromQ(localMakeStarts(nLHS, nJit, d, LBq, UBq, UBq - LBq, q0, lhsOK, tf));
else
    makeStarts = @() localMakeStarts(nLHS, nJit, d, LB, UB, span, z0, lhsOK, false(1, d));
end

starts = makeStarts();
initialguess = starts;
sp = CustomStartPointSet(toQ(starts));

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

while flagg <= 0 && attempt < maxAttempts   % 0 = no local minimum found (a run hit its limits): retry too

    attempt = attempt + 1;

    if attempt > 1
        starts = makeStarts();                 % fresh draws each retry
        initialguess = [initialguess; starts]; %#ok<AGROW>
        sp = CustomStartPointSet(toQ(starts));
    end

    [P, fval, flagg, outpt, allmins] = run(ms, problem, sp);
    P = fromQ(P);

    % Track the best finite local minimum seen across attempts
    if ~isempty(allmins)
        [fv, idx] = min([allmins.Fval]);
        if isfinite(fv) && fv < bestFval
            bestFval = fv;
            bestP = fromQ(allmins(idx).X);
        end
    end

end

if flagg <= 0
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

% Optional derivative-free polish of the returned solution (see header).
% CONVERGED describes the parameter vector actually returned: true when it is the
% MultiStart solution with a positive exit flag and no polish changed it, or when
% the polish pass that produced it converged (exit flag 1) with a feasible raw
% candidate. An improved but unconverged or projected polish point is still
% returned (it is a better point) but is reported as not converged.
multistartExitflag = flagg;
converged = flagg > 0;
polishExitflag = NaN; polishImproved = false; confirmExitflag = NaN;
if isfield(params, 'polish') && ~isempty(params.polish) && params.polish
    freeIdx = find(UBq > LBq);
    if ~isempty(freeIdx) && all(isfinite(P))
        qP = toQ(P(:)');
        nFree = numel(freeIdx);
        optsNM = optimset('TolX', 1e-8, 'TolFun', 1e-8, 'Display', 'off', ...
            'MaxFunEvals', 500 * nFree, 'MaxIter', 500 * nFree);
        for pass = 1:2
            try
                fpol = @(qf) localPolishObjective(qf, qP, freeIdx, LBq, UBq, fq);
                [qfRaw, ~, polishExitflag] = fminsearch(fpol, qP(freeIdx), optsNM);
                qfRaw = qfRaw(:)';
                feasible = all(qfRaw >= LBq(freeIdx)) && all(qfRaw <= UBq(freeIdx));
                qf = min(max(qfRaw, LBq(freeIdx)), UBq(freeIdx));    % project to the box ...
                qTry = qP; qTry(freeIdx) = qf;
                fTry = fq(qTry);                                       % ... and evaluate there before accepting
                improvedNow = isfinite(fTry) && fTry < fval;
                if improvedNow
                    qP = qTry; P = fromQ(qP); fval = fTry; polishImproved = true;
                    converged = (polishExitflag == 1) && feasible;
                elseif polishExitflag == 1 && feasible
                    converged = converged || flagg > 0;                % polish confirms the SQP point
                end
                % A second pass only when the first one left the point unconverged
                % (budget exhausted or projected): it restarts from the improved point.
                if converged || ~improvedNow, break, end
            catch
                break   % keep the current solution if the polish fails for any reason
            end
        end
        % Confirmation step: a polished point that is still unconverged (budget
        % exhausted, or a raw candidate outside the box, which is what happens when
        % the optimum lies ON a bound) is handed to a bounded local optimization
        % from that point; its exit flag then describes the returned vector.
        if polishImproved && ~converged
            try
                optConfirm = optimoptions(options, 'Display', 'off');
                [qc, fc2, flc] = fmincon(fq, qP, [], [], [], [], LBq, UBq, [], optConfirm);
                qc = min(max(qc(:)', LBq), UBq);
                if flc > 0 && isfinite(fc2) && fc2 <= fval + 1e-9 * max(1, abs(fval))
                    qP = qc; P = fromQ(qP); fval = min(fc2, fval); converged = true;
                    confirmExitflag = flc;
                end
            catch
                % keep the polished point, reported as not converged
            end
        end
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

    fitDiagnostics.LB = LB;   % full optimizer box [params.LB, X0 bounds, LBe] (post-flipMask)   % NEW
    fitDiagnostics.UB = UB;   % used downstream for PII_box = log10(UB) - log10(LB)              % NEW
    fitDiagnostics.multistartExitflag = multistartExitflag;
    fitDiagnostics.polishExitflag = polishExitflag;   % NaN = no polish; 1 converged; 0 budget exhausted
    fitDiagnostics.polishImproved = polishImproved;
    fitDiagnostics.confirmExitflag = confirmExitflag; % fmincon exit flag of the confirmation step (NaN = not run)
    fitDiagnostics.converged = converged;             % convergence status of the RETURNED vector
    fitDiagnostics.logscale = useLog;
    fitDiagnostics.scaleProblem = scaleProblem;
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
function starts = localMakeStarts(nLHS, nJit, d, LB, UB, span, z0, lhsOK, logMask)
% logMask(i) = true when coordinate i is a log10 value: its jitter is then additive,
% q0 + N(0, 0.1) (about +/-25% multiplicative in the original scale); the other
% coordinates keep the multiplicative jitter.

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
    if nargin < 9 || isempty(logMask), logMask = false(1, d); end
    jitter = 10.^(0.1*randn(nJit, d));
    S2 = repmat(z0, nJit, 1) .* jitter;
    zeroCols = (z0 == 0) & ~logMask;
    if any(zeroCols)
        S2(:, zeroCols) = repmat(z0(zeroCols), nJit, 1) + ...
            0.01*randn(nJit, sum(zeroCols)) .* repmat(span(zeroCols), nJit, 1);
    end
    if any(logMask)
        S2(:, logMask) = repmat(z0(logMask), nJit, 1) + 0.1*randn(nJit, sum(logMask));
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


% <============================================================================>
% < log10 reparametrization helpers (paramsX.logscale)                        ==>
% <============================================================================>
function q = localToQ(z, tf)
% Rows of z are points in the original coordinates; transformed columns -> log10.
q = z;
if any(tf)
    q(:, tf) = log10(max(z(:, tf), realmin));
end
end

function v = localPolishObjective(qf, qFull, freeIdx, LBq, UBq, fq)
% Objective for the Nelder-Mead polish: free coordinates only, box enforced by penalty.
if any(qf(:)' < LBq(freeIdx)) || any(qf(:)' > UBq(freeIdx))
    v = 1e10;
    return
end
q = qFull;
q(freeIdx) = qf(:)';
v = fq(q);
end

function z = localFromQ(q, tf)
z = q;
if any(tf)
    z(:, tf) = 10 .^ q(:, tf);
end
end
