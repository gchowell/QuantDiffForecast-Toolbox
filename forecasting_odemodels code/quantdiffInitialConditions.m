function IC = quantdiffInitialConditions(z, params, vars)
%QUANTDIFFINITIALCONDITIONS  Initial state vector of the ODE for one parameter vector.
%
%   IC = quantdiffInitialConditions(z, params, vars)
%
%   z       full optimizer vector [params, X0 of the fitted states, alpha, d]
%   params  parameter struct (params.num or params.label gives the number of
%           model parameters)
%   vars    variable struct (vars.initial, vars.fit_index and, optionally,
%           vars.initialFcn and vars.window)
%
%   Without vars.initialFcn (the original behaviour, unchanged):
%
%       IC = vars.initial;  IC(vars.fit_index) = X0;
%
%   With vars.initialFcn the initial conditions may depend on the parameters:
%
%       IC = vars.initialFcn(p, x0, vars);   IC(vars.fit_index) = x0;
%
%   where p = z(1:params.num) are the model parameters of this trial and x0 the
%   initial values of the fitted states (the first observation of each series
%   when params.fixI0 = 1, the trial value when params.fixI0 = 0). The function
%   returns one value per state, in the order of vars.label.
%
%   Why: with vars.initial alone an initial condition cannot follow a parameter.
%   For the SEIR example, S(0) = N - I(0) stayed at its options-file value when
%   N was estimated, so beta and N entered the model only as beta/N and could
%   not be told apart whatever the data. Unobserved initial states (E(0), I(0))
%   could not be estimated either, and later rolling windows restarted them from
%   the same values. Examples:
%
%     % S(0) follows N (states S E I R C; parameters beta kappa gamma N)
%     vars.initialFcn = @(p, x0, vars) [p(4) - 4, 0, 4, 0, x0];
%
%     % I(0) and E(0) as parameters 5 and 6 (the ODE function ignores them)
%     vars.initialFcn = @(p, x0, vars) [p(4)-p(5)-p(6), p(6), p(5), 0, x0];
%
%     % rolling windows: start from a depleted population, counting the cases
%     % reported before the window as removed (an approximation)
%     vars.initialFcn = @(p, x0, vars) [p(4)-p(5)-p(6)-vars.window.priorCum, p(6), p(5), vars.window.priorCum, x0];
%
%   The fitted states are always set to x0 AFTER the function is called. The
%   toolbox relies on that: with params.fixI0 = 1 the first fitted point equals
%   the first observation by construction, which is what quantdiffEffectiveN,
%   the AICc and the bootstrap assume.
%
%   vars.window describes the calibration window and is set by the runners:
%     .tstart    row index of the first calibration row
%     .rows      row indices of the calibration window
%     .time0     time of the first calibration row
%     .firstObs  first observation of every fitted series in the window
%     .priorCum  sum of each observed series over the rows BEFORE the window
%                (zeros for a window that starts at row 1)
%   When it is absent (simulation, model plots) a first-window default is
%   supplied, so a function may always read it.
%
%   An invalid result (wrong length, non-finite, complex or negative values)
%   raises QuantDiffForecast:InvalidInitialFcn. During optimization
%   parameterSearchODE turns that into a failed trial, so a parameter vector
%   whose initial conditions are infeasible (e.g. N below the initial number
%   infected) is simply never selected.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

nFit = numel(vars.fit_index);
if isfield(params, 'num') && ~isempty(params.num)
    nPar = params.num;
else
    nPar = numel(params.label);
end
if numel(z) < nPar + nFit
    error('QuantDiffForecast:InvalidInitialConditions', ...
        'The parameter vector has %d entries; at least %d are needed (%d parameters and %d initial values).', ...
        numel(z), nPar + nFit, nPar, nFit);
end
x0 = reshape(z(nPar + (1:nFit)), 1, nFit);

if isfield(vars, 'initialFcn') && ~isempty(vars.initialFcn)

    if ~isa(vars.initialFcn, 'function_handle')
        error('QuantDiffForecast:InvalidInitialFcn', ...
            'vars.initialFcn must be a function handle, e.g. @(p, x0, vars) [p(4)-4, 0, 4, 0, x0].');
    end

    % First-window default, so the function can always read vars.window.
    if ~isfield(vars, 'window') || isempty(vars.window)
        vars.window = struct('tstart', 1, 'rows', [], 'time0', 0, ...
            'firstObs', x0, 'priorCum', zeros(1, nFit));
    end

    p  = reshape(z(1:nPar), 1, nPar);
    IC = vars.initialFcn(p, x0, vars);

    nStates = numel(vars.initial);
    if ~isnumeric(IC) || ~isvector(IC) || numel(IC) ~= nStates
        error('QuantDiffForecast:InvalidInitialFcn', ...
            'vars.initialFcn must return a numeric vector with %d values, one per state in vars.label.', nStates);
    end
    if ~isreal(IC) || any(~isfinite(IC(:))) || any(IC(:) < 0)
        error('QuantDiffForecast:InvalidInitialFcn', ...
            'vars.initialFcn returned non-finite, complex or negative initial conditions for this parameter vector.');
    end
    IC = reshape(double(IC), size(vars.initial));

else

    IC = vars.initial;

end

% The fitted states always start at x0 (see header).
IC(vars.fit_index) = x0;

end
