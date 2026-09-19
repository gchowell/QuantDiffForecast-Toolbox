function [options, profile] = quantdiffODEOptions(IC)
%QUANTDIFFODEOPTIONS One numerical policy for estimation and reconstruction.
% Preserve the original objective's speed-oriented tolerances. All fitting,
% final-fit reconstruction and bootstrap forecasting use this helper.
% Independent synthetic-data generation in Run_simulate_ODEModel is separate.
% MaxStep is deliberately unspecified: MATLAB chooses its default.
% As in the original toolbox, all state variables are constrained nonnegative.
% Models with signed states require an explicit, consistently applied policy.

if ~isnumeric(IC) || ~isvector(IC) || isempty(IC) || ...
        ~isreal(IC) || any(~isfinite(IC(:)))
    error('QuantDiffForecast:InvalidInitialConditions', ...
        'Initial conditions must be a finite real numeric vector.');
end
options = odeset('RelTol',1e-6, 'AbsTol',1e-8, ...
    'NonNegative',1:numel(IC));
if nargout > 1
    profile = struct('id','quantdiff-consistent-v1', ...
        'solver','ode15s', 'RelTol',options.RelTol, ...
        'AbsTol',options.AbsTol, 'MaxStep',options.MaxStep, ...
        'NonNegative',options.NonNegative);
end
end
