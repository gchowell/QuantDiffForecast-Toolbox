"""Observation noise: AddErrorStructure.m, estimatePearsonDispersion.m,
estimateLaplaceScale.m and quantdiffNegativeNoiseShare.m.

MATLAB's normrnd/poissrnd/nbinrnd/rand used the global random stream; here
every function takes a numpy ``Generator`` explicitly. numpy's
``negative_binomial(n, p)`` has the same (r, p) parameterization as
MATLAB's ``nbinrnd(r, p)`` (failures before the r-th success).
"""
from __future__ import annotations

import math

import numpy as np


def add_error_structure(yi, M, dist1, factor1, d, rng):
    """Draw M noisy observation curves around a fitted curve.

    yi is the *cumulative* curve (the runners pass cumsum(fitcurve)), so the
    increments diff(yi) are the fitted observations. Row 0 of every draw is
    yi[0] (the first observation is never perturbed). Returns len(yi) x M.
    """
    yi = np.asarray(yi, dtype=float).ravel()
    n = yi.size
    curves = np.zeros((n, M))
    curves[0, :] = yi[0]
    if n < 2 or M == 0:
        return curves
    inc = np.diff(yi)[:, None] * np.ones((1, M))
    lam = np.abs(inc)
    eps1 = 0.001
    if dist1 == 0:
        draws = rng.normal(inc, factor1)
    elif dist1 == 1:
        draws = rng.poisson(lam)
    elif dist1 in (2, 3, 4, 5):
        mean1 = np.where(lam == 0, eps1, lam)
        if dist1 == 2:
            var1 = mean1 * factor1
        elif dist1 == 3:
            var1 = mean1 + mean1 * factor1
        elif dist1 == 4:
            var1 = mean1 + factor1 * mean1 ** 2
        else:
            var1 = mean1 + factor1 * mean1 ** d
        p1 = mean1 / var1
        r1 = mean1 * p1 / (1 - p1)
        draws = rng.negative_binomial(r1, p1).astype(float)
    elif dist1 == 6:
        U = rng.random(lam.shape) - 0.5
        draws = lam - factor1 * np.sign(U) * np.log(1 - 2 * np.abs(U))
    else:
        raise ValueError(f"Unsupported error structure dist1 = {dist1}")
    curves[1:, :] = draws
    return curves


def estimate_pearson_dispersion(ydata, fitcurve, numModelParams, nSeries=1, fixI0=0, minMean=0.5):
    y = np.asarray(ydata, dtype=float).ravel()
    mu = np.asarray(fitcurve, dtype=float).ravel()
    nT = y.size // nSeries
    use = np.isfinite(y) & np.isfinite(mu)
    fixed = np.zeros(y.size, dtype=bool)
    if fixI0:
        fixed[::nT] = True
    small = mu < minMean
    use = use & ~fixed & ~small
    dof = int(use.sum()) - numModelParams
    phi = np.nan if dof < 1 else float(np.sum((y[use] - mu[use]) ** 2 / mu[use]) / dof)
    info = {"nUsed": int(use.sum()), "dof": dof, "minMean": minMean,
            "excludedSmallMean": int((small & ~fixed).sum()), "excludedFixed": int(fixed.sum())}
    return phi, info


def estimate_laplace_scale(ydata, fitcurve, numModelParams, nSeries=1, fixI0=0):
    y = np.asarray(ydata, dtype=float).ravel()
    mu = np.asarray(fitcurve, dtype=float).ravel()
    nT = y.size // nSeries
    use = np.isfinite(y) & np.isfinite(mu)
    fixed = np.zeros(y.size, dtype=bool)
    if fixI0:
        fixed[::nT] = True
    use = use & ~fixed
    dof = use.sum() - numModelParams / 2
    b = np.nan if use.sum() <= numModelParams else float(np.sum(np.abs(y[use] - mu[use])) / dof)
    return b, {"nUsed": int(use.sum()), "dof": dof, "excludedFixed": int(fixed.sum())}


def negative_noise_share(fitcurve, scale, distBoot, nSeries=1):
    """Expected share of negative bootstrap draws for normal / Laplace noise."""
    info = {"nPoints": 0, "minMean": np.nan, "worstPoint": np.nan}
    if distBoot not in (0, 6) or not np.isscalar(scale) or not np.isfinite(scale) or scale <= 0:
        return np.nan, info
    mu = np.abs(np.asarray(fitcurve, dtype=float).ravel())
    nT = mu.size // nSeries
    noisy = np.ones(mu.size, dtype=bool)
    noisy[::nT] = False
    noisy &= np.isfinite(mu)
    if not noisy.any():
        return np.nan, info
    mu = mu[noisy]
    if distBoot == 0:
        p = 0.5 * np.array([math.erfc(m / (scale * math.sqrt(2))) for m in mu])
    else:
        p = 0.5 * np.exp(-mu / scale)
    info = {"nPoints": mu.size, "minMean": float(mu.min()), "worstPoint": float(p.max())}
    return float(p.mean()), info
