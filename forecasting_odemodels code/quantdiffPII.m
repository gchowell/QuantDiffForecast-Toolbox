function [PII, rangeFrac, scale, PIIbox, ratio] = quantdiffPII(lo, hi, LBbox, UBbox, scaleUser)
%QUANTDIFFPII  Parameter identifiability index and its companions, defined for every parameter.
%
%   [PII, rangeFrac, scale, PIIbox, ratio] = quantdiffPII(lo, hi, LBbox, UBbox)
%   [...] = quantdiffPII(lo, hi, LBbox, UBbox, scaleUser)
%
%   lo, hi        2.5% and 97.5% bootstrap quantiles of each parameter (row vectors)
%   LBbox, UBbox  optimizer bounds of the same parameters (fitDiagnostics.LB / .UB)
%   scaleUser     optional override of SCALE per parameter (params.PIIscale):
%                 1 = log scale, 0 = range fraction, NaN = automatic. Shorter
%                 vectors are padded with NaN.
%
%   PII        log10(hi/lo), the width of the 95% interval in decades (1 = the
%              interval spans a factor of 10). Defined only when lo > 0;
%              NaN otherwise. Formerly a zero lower quantile gave Inf and a
%              negative one a complex number, and both entered medians and plots.
%   rangeFrac  (hi - lo)/(UBbox - LBbox): the share of the search range that the
%              interval still covers, in [0,1]. Always defined except for pinned
%              parameters (UBbox == LBbox -> NaN). 1 = the data add nothing to
%              the bounds. It depends on how wide the bounds were set, so read
%              it together with PII, not instead of it.
%   scale      which measure is the headline for each parameter:
%                1  log scale (PII): the lower bound is positive, so relative
%                   precision is meaningful (rates, population sizes, alpha)
%                0  range fraction: the lower bound is <= 0 (proportions,
%                   anything that may be zero), where a log ratio is undefined
%                   or inflated
%              This is the rule params.logscale uses for the search scale.
%   PIIbox     log10(UBbox/LBbox), the largest PII the bounds allow; NaN when
%              LBbox <= 0 (formerly Inf, which made the ratio a false zero).
%   ratio      bound-censoring diagnostic, ~1 = the interval is the search range:
%                scale 1:  PII ./ PIIbox   (NaN for pinned parameters)
%                scale 0:  rangeFrac
%
%   PII is unchanged for every parameter whose lower quantile and lower bound
%   are positive.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

lo = double(lo(:)'); hi = double(hi(:)'); LBbox = double(LBbox(:)'); UBbox = double(UBbox(:)');
n = numel(lo);
if numel(hi) ~= n || numel(LBbox) ~= n || numel(UBbox) ~= n
    error('quantdiffPII:size', 'lo, hi, LBbox and UBbox must have the same number of elements.');
end
if ~isreal(lo) || ~isreal(hi) || ~isreal(LBbox) || ~isreal(UBbox)
    error('quantdiffPII:complex', 'Inputs must be real.');
end

% Log-scale index: only where both quantiles are positive.
PII = nan(1, n);
ok = lo > 0 & hi > 0;
PII(ok) = log10(hi(ok)) - log10(lo(ok));

% Share of the search range covered by the interval.
width = UBbox - LBbox;
rangeFrac = nan(1, n);
free = width > 0 & isfinite(width) & isfinite(lo) & isfinite(hi);
rangeFrac(free) = min(max((hi(free) - lo(free)) ./ width(free), 0), 1);

% Headline scale: automatic rule, then the user's override where given.
scale = double(LBbox > 0);
if nargin >= 5 && ~isempty(scaleUser)
    scaleUser = double(scaleUser(:)');
    m = min(numel(scaleUser), n);
    given = false(1, n);
    given(1:m) = isfinite(scaleUser(1:m));
    if any(scaleUser(given(1:m)) ~= 0 & scaleUser(given(1:m)) ~= 1)
        error('quantdiffPII:scale', 'params.PIIscale entries must be 1 (log scale), 0 (range fraction) or NaN (automatic).');
    end
    scale(given) = scaleUser(given);
end

% Largest PII the bounds allow; undefined for a nonpositive lower bound.
PIIbox = nan(1, n);
pos = LBbox > 0 & UBbox > 0;
PIIbox(pos) = log10(UBbox(pos)) - log10(LBbox(pos));

% Bound-censoring diagnostic on the headline scale of each parameter.
ratio = nan(1, n);
isLog = scale == 1;
ratio(isLog)  = PII(isLog) ./ PIIbox(isLog);
ratio(~isLog) = rangeFrac(~isLog);
ratio(isLog & PIIbox == 0) = NaN;          % pinned parameter: 0/0, made explicit
end
