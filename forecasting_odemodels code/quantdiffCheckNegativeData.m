function quantdiffCheckNegativeData(method1, data, rows, caller)
%QUANTDIFFCHECKNEGATIVEDATA  Stop when a count likelihood is asked to fit negative observations.
%
%   quantdiffCheckNegativeData(method1, data, rows)
%   quantdiffCheckNegativeData(method1, data, rows, caller)
%
%   method1  estimation method (global method1 of the options file)
%   data     the data matrix as read from the input file: column 1 = time,
%            columns 2..end = one observed series each
%   rows     row numbers of DATA used for calibration (all rolling windows:
%            tstart1 : tend1 + windowsize1 - 1). Rows used only to evaluate
%            forecasts are not checked.
%   caller   name shown in the message (default 'QuantDiffForecast')
%
%   The Poisson and negative binomial likelihoods (method1 = 1, 3, 4, 5) give
%   probability zero to a negative count, so a negative observation has no
%   likelihood. Before this check the objective still returned a number that
%   treated the same point in two different ways (clipped to zero in the gamma
%   terms, used as is in the y*log terms), and the fit, AICc and profile
%   intervals inherited it. The run now stops before any fitting, lists the
%   offending rows and says what to do.
%
%   Least squares (method1 = 0) and the sum of absolute deviations
%   (method1 = 6) accept negative values, so nothing is checked for them.
%
%   Negative values usually come from daily increments of a cumulative series
%   that was revised downward. Correct the data first: set the value to zero,
%   or better, spread the correction over the earlier dates it belongs to.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 4 || isempty(caller), caller = 'QuantDiffForecast'; end

if ~any(method1 == [1 3 4 5])
    return                                  % least squares / SAD: negative values are allowed
end

rows = rows(:);
rows = rows(rows >= 1 & rows <= size(data, 1) & rows == fix(rows));   % the runners validate the window range themselves
vals = data(rows, 2:end);
[rr, cc] = find(vals < 0);
if isempty(rr)
    return
end

% List the offending points in row order (at most maxShow of them).
[~, order] = sortrows([rows(rr) cc]);
rr = rr(order); cc = cc(order);
nBad = numel(rr);
maxShow = 10;
nShow = min(nBad, maxShow);
lines = cell(1, nShow);
for k = 1:nShow
    lines{k} = sprintf('  row %d (time %g), series %d: %g', ...
        rows(rr(k)), data(rows(rr(k)), 1), cc(k), vals(rr(k), cc(k)));
end
more = '';
if nBad > nShow
    more = sprintf('\n  ... and %d more', nBad - nShow);
end

switch method1
    case 1
        mname = 'Poisson';
    otherwise
        mname = 'negative binomial';
end

error('QuantDiffForecast:NegativeCounts', ...
    ['%s: the calibration data contain %d negative value(s), which the %s likelihood ' ...
     '(method1 = %d) cannot fit: a negative count has probability zero.\n%s%s\n' ...
     'Correct the data before fitting: set these values to zero, or spread the correction ' ...
     'over the earlier dates it belongs to (negative daily counts usually come from a ' ...
     'cumulative series that was revised downward). Least squares (method1 = 0) and SAD ' ...
     '(method1 = 6) accept negative values.'], ...
    caller, nBad, mname, method1, strjoin(lines, sprintf('\n')), more);
end
