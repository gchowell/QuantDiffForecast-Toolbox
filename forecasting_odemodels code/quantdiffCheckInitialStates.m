function flagged = quantdiffCheckInitialStates(vars, lastWindowStart, caller)
%QUANTDIFFCHECKINITIALSTATES  Warn when a later window would restart the unobserved states.
%
%   flagged = quantdiffCheckInitialStates(vars, lastWindowStart)
%   flagged = quantdiffCheckInitialStates(vars, lastWindowStart, caller)
%
%   vars             variable struct of the options file
%   lastWindowStart  row index of the last calibration-window start (tend1 in
%                    the runners, tstart1 in the profile likelihood)
%   caller           name shown in the message
%
%   A calibration window that starts after row 1 sets the fitted states to the
%   window's first observation, but every other state starts from vars.initial,
%   the values written for the START of the series. For a compartmental model
%   that means each later window restarts the epidemic (for the SEIR example:
%   I(0) = 4 and a fully susceptible population) however many cases came before,
%   and the parameter estimates absorb the mismatch.
%
%   The warning is issued when lastWindowStart > 1, the model has states that are
%   not fitted, and the options file has no vars.initialFcn. Supplying
%   vars.initialFcn (which can read vars.window.priorCum; see
%   quantdiffInitialConditions) or using expanding windows that all start at
%   row 1 removes it. Models whose only state is the fitted one (EXP, GGM, GLM)
%   are not flagged.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 3 || isempty(caller), caller = 'QuantDiffForecast'; end

hasFcn  = isfield(vars, 'initialFcn') && ~isempty(vars.initialFcn);
nStates = numel(vars.initial);
flagged = lastWindowStart > 1 && nStates > numel(vars.fit_index) && ~hasFcn;

if flagged
    warning('QuantDiffForecast:rollingInitialStates', ...
        ['%s: calibration windows start as late as row %d, but the %d state(s) that are not fitted ' ...
         'start every window from vars.initial, the values given for the start of the series. ' ...
         'Set vars.initialFcn so the initial conditions follow the window (vars.window.priorCum holds ' ...
         'the counts observed before it; see quantdiffInitialConditions), or use windows that all ' ...
         'start at row 1. Silence with warning(''off'',''QuantDiffForecast:rollingInitialStates'').'], ...
        caller, lastWindowStart, nStates - numel(vars.fit_index));
end
end
