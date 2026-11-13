function quantdiffCheckMethod(method1, caller)
%QUANTDIFFCHECKMETHOD  Stop with a clear message when method1 is not a supported estimation method.
%
%   quantdiffCheckMethod(method1)
%   quantdiffCheckMethod(method1, caller)     caller = name shown in the message
%
%   Supported estimation methods (global method1, set in the options file):
%     0  nonlinear least squares
%     1  maximum likelihood, Poisson
%     3  maximum likelihood, negative binomial, VAR = mean + alpha*mean
%     4  maximum likelihood, negative binomial, VAR = mean + alpha*mean^2
%     5  maximum likelihood, negative binomial, VAR = mean + alpha*mean^d
%     6  sum of absolute deviations (Laplace)
%
%   method1 = 2 is NOT an estimation method. (2 is only a value of dist1: the
%   negative binomial bootstrap noise VAR = factor*mean used with method1 = 0.)
%   Without this check an unsupported value failed later inside fit_model with
%   "Undefined function or variable 'LBe'".
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 2 || isempty(caller), caller = 'QuantDiffForecast'; end
supported = [0 1 3 4 5 6];

if isempty(method1)
    error('QuantDiffForecast:UnsupportedMethod', ...
        ['%s: method1 is empty. The options file must declare "global method1" and set it to ' ...
         'one of 0, 1, 3, 4, 5 or 6.'], caller);
end
if ~isnumeric(method1) || ~isscalar(method1) || ~isreal(method1) || ~any(method1 == supported)
    if isnumeric(method1) && isscalar(method1) && isreal(method1)
        shown = sprintf('%g', method1);
    else
        shown = sprintf('a %s of size %s', class(method1), mat2str(size(method1)));
    end
    hint = '';
    if isnumeric(method1) && isscalar(method1) && isreal(method1) && method1 == 2
        hint = [' For least squares with negative binomial bootstrap noise (VAR = factor*mean) ' ...
                'set method1 = 0 and dist1 = 2.'];
    end
    error('QuantDiffForecast:UnsupportedMethod', ...
        ['%s: method1 = %s is not supported; use 0 (least squares), 1 (Poisson), ' ...
         '3, 4 or 5 (negative binomial) or 6 (SAD/Laplace).%s'], caller, shown, hint);
end
end
