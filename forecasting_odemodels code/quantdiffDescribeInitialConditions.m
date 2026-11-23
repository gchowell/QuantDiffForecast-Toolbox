function txt = quantdiffDescribeInitialConditions(params, vars)
%QUANTDIFFDESCRIBEINITIALCONDITIONS  One line for the settings summaries: the initial conditions in use.
%
%   txt = quantdiffDescribeInitialConditions(params, vars)
%
%   Without vars.initialFcn the line shows vars.initial, as before:
%
%       Initial Conditions: 549996  0  4  0  4
%
%   With vars.initialFcn it shows the values the function gives at
%   params.initial and says where they come from, because vars.initial is then
%   only the fallback and the conditions change with the parameters:
%
%       Initial Conditions (vars.initialFcn at params.initial; they follow the parameters): ...
%
%   Used by the summaries printed by plotODEModel, plotFit_ODEModel and
%   Run_simulate_ODEModel.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if isfield(vars, 'initialFcn') && ~isempty(vars.initialFcn)
    z  = [params.initial(:)' reshape(vars.initial(vars.fit_index), 1, []) 1 1];
    IC = quantdiffInitialConditions(z, params, vars);
    txt = ['Initial Conditions (vars.initialFcn at params.initial; they follow the parameters): ', num2str(IC(:)')];
else
    txt = ['Initial Conditions: ', num2str(vars.initial)];
end
end
