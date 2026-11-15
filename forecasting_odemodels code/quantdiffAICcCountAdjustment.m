function [adj, info] = quantdiffAICcCountAdjustment(method1, ydata, nSeries, fixI0, alpha, d)
%QUANTDIFFAICCCOUNTADJUSTMENT  Turn the count objective into a complete log-likelihood for AICc.
%
%   [adj, info] = quantdiffAICcCountAdjustment(method1, ydata, nSeries, fixI0, alpha, d)
%
%   For the Poisson and negative binomial methods (method1 = 1, 3, 4, 5) the
%   fitting objective fval is the negative log-likelihood WITHOUT the term that
%   depends only on the data, sum(log(y!)). That is harmless for fitting, the
%   bootstrap and the profile likelihood (it cancels in differences), but it
%   made the AICc of these methods incomparable with the AICc of least squares
%   (method1 = 0) and SAD (method1 = 6), whose formulas in getAICc are complete
%   Gaussian / Laplace log-likelihoods: for an outbreak of a few thousand cases
%   the missing term is tens of thousands of AICc units.
%
%   ADJ is the amount to add to fval so that fval + ADJ is the full negative
%   log-likelihood of the observations that count toward the sample size
%   (quantdiffEffectiveN):
%
%       adj = sum(gammaln(y + 1)) over the informative observations
%             - (objective of the fixed first points)          [fixI0 = 1 only]
%
%   With params.fixI0 = 1 the first point of every series is matched exactly by
%   construction; it is not counted in n, so its likelihood term is taken out
%   as well (getAICc then uses the same observations in the likelihood and in
%   the small-sample correction, as it already does for methods 0 and 6). The
%   term is evaluated with quantdiffObjectiveValue at mean = observation, with
%   the fitted alpha and d, i.e. exactly what fval contains for those points.
%
%   ADJ is 0 for methods 0 and 6.
%
%   Inputs
%     method1   estimation method
%     ydata     stacked observations (one block of equal length per series)
%     nSeries   number of fitted series, length(vars.fit_index)
%     fixI0     params.fixI0
%     alpha, d  fitted error-model parameters (ignored for method1 = 1)
%
%   info: struct with logFactorialTerm, firstPointTerm, nInformative.
%
%   fval itself is not changed anywhere: only the AICc value moves. For a given
%   window the shift is the same for every Poisson model, so Poisson model
%   rankings are unchanged; for negative binomial models the first-point term
%   depends on each model's alpha, so differences can move slightly.
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

adj = 0;
info = struct('logFactorialTerm', 0, 'firstPointTerm', 0, 'nInformative', numel(ydata));
if ~any(method1 == [1 3 4 5])
    return
end
if nargin < 5 || isempty(alpha), alpha = 1; end
if nargin < 6 || isempty(d),     d = 1;     end

y = double(ydata(:));
if mod(numel(y), nSeries) ~= 0
    error('quantdiffAICcCountAdjustment:series', ...
        '%d observations cannot be split into %d series of equal length.', numel(y), nSeries);
end
nT = numel(y) / nSeries;

pinned = false(size(y));
if ~isempty(fixI0) && fixI0 == 1
    pinned(1:nT:end) = true;                 % first point of every series
end

yInf = max(y(~pinned), 0);                   % negative counts are refused upstream (quantdiffCheckNegativeData)
logFact = sum(gammaln(yInf + 1));

firstTerm = 0;
if any(pinned)
    y1 = y(pinned);
    mu1 = y1;
    mu1(mu1 == 0) = 0.001;                   % the objective's own replacement for a zero mean
    firstTerm = quantdiffObjectiveValue(mu1, y1, method1, alpha, d);
end

adj = logFact - firstTerm;
info = struct('logFactorialTerm', logFact, 'firstPointTerm', firstTerm, 'nInformative', nnz(~pinned));
end
