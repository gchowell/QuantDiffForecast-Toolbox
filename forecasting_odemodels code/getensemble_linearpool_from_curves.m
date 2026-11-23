function [Qens, forecast1, Qmodels, meta, curvesEns] = getensemble_linearpool_from_curves(weights1, varargin)
% Weighted Predictive Distribution Ensemble (Exact Empirical Linear Pool).
%   - Inputs: K model matrices A1..AK, each T x M (same T, M), plus weights1.
%   - Internally computes per-model empirical quantiles at specified q-levels.
%   - Forms F_ens(x) = sum_i w_i F_i(x) exactly from the empirical draws and
%     evaluates the generalized inverse Q(p) = inf{x:F_ens(x) >= p}.
%
% Usage:
%   [Qens, forecast1] = getensemble_linearpool_from_curves(w, A1, A2, A3);
%   [Qens, forecast1] = getensemble_linearpool_from_curves(w, A1, A2, ...
%       'QLevels', 0.01:0.01:0.99, ...
%       'ProbsOut', [0.5 0.025 0.975]);
%
% Outputs:
%   Qens       : T x Lout matrix of exact ensemble quantiles at 'ProbsOut'.
%   forecast1  : T x 3 matrix [median, LB_2.5%, UB_97.5%].
%   Qmodels    : 1 x K cell, each T x L matrix of model empirical quantiles.
%   meta       : struct with fields including QLevels, ProbsOut, K, T, M,
%                Tol, MaxIter, weights, ClampNonneg, and method details.
%   curvesEns  : Optional T x Nens reproducible sample of complete trajectories
%                from the model mixture. Within-trajectory temporal dependence
%                is preserved. If only four outputs are requested, trajectory
%                sampling is skipped.
%
% Name-value options retained for compatibility:
%   'QLevels'      Default 0.01:0.01:0.99. Strictly increasing values in [0,1].
%   'ProbsOut'     Default [0.5 0.025 0.975]. Values in [0,1].
%   'Tol'          Default 1e-6. Retained for old calls but not used because
%                  the exact empirical calculation requires no inversion.
%   'MaxIter'      Default 60. Retained for old calls but not used.
%   'ClampNonneg'  Default true. If true, every draw is clamped before any
%                  component or ensemble quantile is calculated.
%   'Nens'         Optional number of trajectories in curvesEns. Default is M,
%                  the number of trajectories in each input model. Set to zero
%                  to return an empty curvesEns.
%   'Seed'         Nonnegative integer seed for curvesEns. Default 1. A local
%                  random stream is used, so the global RNG state is unchanged.
%
% Important QuantDiffForecast note:
%   Qens contains pointwise quantiles, not predictive draws. Use curvesEns,
%   rather than Qens, with a routine that expects bootstrap/predictive draws.
%   Call this function separately for process curves and for observation-error-
%   augmented predictive curves. Use predictive curves for WIS, coverage, and
%   prediction intervals. Request curvesEns when complete temporal trajectories
%   or path-dependent outcomes are required.

    % ---------- split varargin into model matrices vs name-value ----------
    isName = cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), varargin);
    firstNV = find(isName, 1, 'first');
    if isempty(firstNV)
        mats = varargin;
        NV = {};
    else
        mats = varargin(1:firstNV-1);
        NV = varargin(firstNV:end);
    end

    % ------------------------ validate weights ---------------------------
    validateattributes(weights1, {'numeric'}, ...
        {'vector','real','finite','nonnegative','nonempty'}, ...
        mfilename, 'weights1', 1);

    weightsInput = double(weights1(:));
    K = numel(mats);
    if K < 1
        error([mfilename ':NoModels'], ...
            'Provide at least one model matrix (T x M).');
    end
    if numel(weightsInput) ~= K
        error([mfilename ':WeightModelMismatch'], ...
            ['weights1 length (%d) must equal the number of model ' ...
             'matrices passed (%d).'], numel(weightsInput), K);
    end

    weightSum = sum(weightsInput);
    if weightSum <= 0 || abs(weightSum - 1) > 1e-8
        error([mfilename ':InvalidWeightSum'], ...
            'weights1 must sum to one (tolerance 1e-8).');
    end
    weightsNormalized = weightsInput / weightSum;
    activeModels = find(weightsNormalized > 0);

    % ------------------------ parse name-value options -------------------
    p = inputParser;
    p.FunctionName = mfilename;
    p.CaseSensitive = false;
    p.PartialMatching = false;
    p.KeepUnmatched = false;

    addParameter(p, 'QLevels', 0.01:0.01:0.99, @is_probability_grid);
    addParameter(p, 'ProbsOut', [0.5 0.025 0.975], @is_probability_vector);
    addParameter(p, 'Tol', 1e-6, @is_positive_finite_scalar);
    addParameter(p, 'MaxIter', 60, @is_valid_iteration_count);
    addParameter(p, 'ClampNonneg', true, @is_scalar_logical_value);
    addParameter(p, 'Nens', [], @is_optional_nonnegative_integer);
    addParameter(p, 'Seed', 1, @is_valid_seed);
    parse(p, NV{:});

    qlevels = double(p.Results.QLevels(:).');
    ProbsOut = double(p.Results.ProbsOut(:).');
    Tol = double(p.Results.Tol);
    MaxIter = double(p.Results.MaxIter);
    ClampNonneg = logical(p.Results.ClampNonneg);
    Nens = p.Results.Nens;
    Seed = double(p.Results.Seed);

    % --------------------- validate and prepare matrices -----------------
    for i = 1:K
        validateattributes(mats{i}, {'numeric'}, ...
            {'2d','real','finite','nonempty'}, mfilename, ...
            sprintf('model matrix %d', i), i + 1);
    end

    [T, M] = size(mats{1});
    if T < 1 || M < 1
        error([mfilename ':EmptyModelDimension'], ...
            'Model matrices must contain at least one time and one draw.');
    end

    for i = 1:K
        [Ti, Mi] = size(mats{i});
        if Ti ~= T || Mi ~= M
            error([mfilename ':ModelSizeMismatch'], ...
                'All model matrices must have the same size T x M.');
        end
        if ClampNonneg
            mats{i} = max(mats{i}, 0);
        end
    end

    if isempty(Nens)
        Nens = M;
    else
        Nens = double(Nens);
    end

    % ----------- per-model empirical quantiles for diagnostics -----------
    Qmodels = cell(1, K);
    equalDrawWeights = ones(1, M) / M;
    for i = 1:K
        Qi = zeros(T, numel(qlevels));
        for t = 1:T
            Qi(t,:) = weighted_empirical_quantiles( ...
                mats{i}(t,:), equalDrawWeights, qlevels);
        end
        Qmodels{i} = Qi;
    end

    % ------ exact weighted empirical linear pool at each forecast time ----
    summaryProbs = [0.5 0.025 0.975];
    requestedTogether = [ProbsOut, summaryProbs];
    [allProbs, ~, probabilityMap] = unique(requestedTogether, 'sorted');
    mapQens = probabilityMap(1:numel(ProbsOut));
    mapSummary = probabilityMap(numel(ProbsOut)+1:end);

    pooledWeights = cell(1, numel(activeModels));
    for a = 1:numel(activeModels)
        i = activeModels(a);
        pooledWeights{a} = repmat(weightsNormalized(i) / M, 1, M);
    end
    pooledWeights = [pooledWeights{:}];
    pooledWeights = pooledWeights / sum(pooledWeights);

    Qall = zeros(T, numel(allProbs));
    for t = 1:T
        pooledValues = cell(1, numel(activeModels));
        for a = 1:numel(activeModels)
            i = activeModels(a);
            pooledValues{a} = mats{i}(t,:);
        end
        pooledValues = [pooledValues{:}];
        Qall(t,:) = weighted_empirical_quantiles( ...
            pooledValues, pooledWeights, allProbs);
    end

    Qens = Qall(:, mapQens);
    forecast1 = Qall(:, mapSummary);

    % -------- optional reproducible whole-trajectory mixture sample -----
    curvesEns = [];
    modelLabels = [];
    selectedIndices = [];
    realizedWeights = [];

    if nargout >= 5 && Nens > 0
        stream = RandStream('mt19937ar', 'Seed', Seed);
        activeWeights = weightsNormalized(activeModels);
        activeWeights = activeWeights / sum(activeWeights);
        cumulativeModelWeights = cumsum(activeWeights);
        cumulativeModelWeights(end) = 1;

        u = rand(stream, 1, Nens);
        if numel(activeModels) == 1
            activeLabels = ones(1, Nens);
        else
            activeLabels = 1 + sum(bsxfun(@gt, u, ...
                cumulativeModelWeights(1:end-1)), 1);
        end

        curvesEns = zeros(T, Nens, 'like', mats{activeModels(1)});
        selectedIndices = zeros(1, Nens);

        for a = 1:numel(activeModels)
            locations = find(activeLabels == a);
            if isempty(locations)
                continue;
            end
            chosen = randi(stream, M, 1, numel(locations));
            modelIndex = activeModels(a);
            curvesEns(:, locations) = mats{modelIndex}(:, chosen);
            selectedIndices(locations) = chosen;
        end

        modelLabels = activeModels(activeLabels).';
        counts = accumarray(modelLabels(:), 1, [K 1]);
        realizedWeights = counts / Nens;
    end

    % ----------------------------- metadata -----------------------------
    meta = struct();
    meta.QLevels = qlevels;
    meta.ProbsOut = ProbsOut;
    meta.K = K;
    meta.T = T;
    meta.M = M;
    meta.Tol = Tol;
    meta.MaxIter = MaxIter;
    meta.TolUsed = false;
    meta.MaxIterUsed = false;
    meta.ClampNonneg = ClampNonneg;
    meta.weightsInput = weightsInput;
    meta.weightsNormalized = weightsNormalized;
    meta.activeModels = activeModels;
    meta.method = 'Exact weighted empirical linear probability pool';
    meta.quantileDefinition = 'inf{x:F(x)>=p}; p=0 min; p=1 max';
    meta.QLevelsUsedForPooling = false;
    meta.QensContainsQuantilesNotDraws = true;
    meta.Nens = Nens;
    meta.Seed = Seed;
    meta.trajectoryEnsembleReturned = ~isempty(curvesEns);
    meta.modelLabels = modelLabels;
    meta.selectedIndices = selectedIndices;
    meta.weightsRealizedInTrajectorySample = realizedWeights;
    meta.curvesEnsContainsWholeTrajectories = true;
    meta.inputRole = ['Not encoded: call separately with process curves and ' ...
                      'predictive curves.'];
end

% ============================== helpers ================================

function q = weighted_empirical_quantiles(values, sampleWeights, probs)
% Generalized-inverse quantiles of a weighted empirical distribution.

    values = double(values(:));
    sampleWeights = double(sampleWeights(:));

    [sortedValues, order] = sort(values, 'ascend');
    sortedWeights = sampleWeights(order);
    sortedWeights = sortedWeights / sum(sortedWeights);
    cumulativeWeights = cumsum(sortedWeights);
    cumulativeWeights(end) = 1; % protect p=1 from floating-point drift

    q = zeros(1, numel(probs));
    for j = 1:numel(probs)
        probability = probs(j);
        if probability <= 0
            index = 1;
        elseif probability >= 1
            index = numel(sortedValues);
        else
            index = find(cumulativeWeights >= probability, 1, 'first');
        end
        q(j) = sortedValues(index);
    end
end

function tf = is_probability_grid(x)
    tf = isnumeric(x) && isvector(x) && ~isempty(x) && isreal(x) && ...
         all(isfinite(x(:))) && all(x(:) >= 0) && all(x(:) <= 1) && ...
         all(diff(double(x(:))) > 0);
end

function tf = is_probability_vector(x)
    tf = isnumeric(x) && isvector(x) && ~isempty(x) && isreal(x) && ...
         all(isfinite(x(:))) && all(x(:) >= 0) && all(x(:) <= 1);
end

function tf = is_positive_finite_scalar(x)
    tf = isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0;
end

function tf = is_valid_iteration_count(x)
    tf = isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && ...
         x >= 10 && x == floor(x);
end

function tf = is_scalar_logical_value(x)
    tf = isscalar(x) && (islogical(x) || ...
        (isnumeric(x) && isreal(x) && isfinite(x) && ismember(x, [0 1])));
end

function tf = is_optional_nonnegative_integer(x)
    tf = isempty(x) || (isnumeric(x) && isscalar(x) && isreal(x) && ...
         isfinite(x) && x >= 0 && x == floor(x));
end

function tf = is_valid_seed(x)
    tf = isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && ...
         x >= 0 && x <= 2^32 - 1 && x == floor(x);
end