"""Interval-width and bound diagnostics: quantdiffPII.m, boundHitFractions.m
and the quantileMCSE local function of the runners."""
from __future__ import annotations

import numpy as np

from .scoring import mquantile


def pii(lo, hi, LBbox, UBbox, scaleUser=None):
    """Returns (PII, rangeFrac, scale, PIIbox, ratio), one value per slot."""
    lo, hi = np.asarray(lo, float).ravel(), np.asarray(hi, float).ravel()
    LBbox, UBbox = np.asarray(LBbox, float).ravel(), np.asarray(UBbox, float).ravel()
    n = lo.size
    PII = np.full(n, np.nan)
    ok = (lo > 0) & (hi > 0)
    PII[ok] = np.log10(hi[ok]) - np.log10(lo[ok])
    width = UBbox - LBbox
    rangeFrac = np.full(n, np.nan)
    free = (width > 0) & np.isfinite(width) & np.isfinite(lo) & np.isfinite(hi)
    rangeFrac[free] = np.clip((hi[free] - lo[free]) / width[free], 0, 1)
    scale = (LBbox > 0).astype(float)
    if scaleUser is not None:
        su = np.asarray(scaleUser, float).ravel()
        m = min(su.size, n)
        given = np.zeros(n, bool)
        given[:m] = np.isfinite(su[:m])
        if np.any((su[:m][given[:m]] != 0) & (su[:m][given[:m]] != 1)):
            raise ValueError("params.PIIscale entries must be 1, 0 or NaN.")
        scale[given] = su[:m][given[:m]]
    PIIbox = np.full(n, np.nan)
    pos = (LBbox > 0) & (UBbox > 0)
    PIIbox[pos] = np.log10(UBbox[pos]) - np.log10(LBbox[pos])
    ratio = np.full(n, np.nan)
    isLog = scale == 1
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio[isLog] = PII[isLog] / PIIbox[isLog]
    ratio[~isLog] = rangeFrac[~isLog]
    ratio[isLog & (PIIbox == 0)] = np.nan
    return PII, rangeFrac, scale, PIIbox, ratio


def bound_hit_fractions(Phatss, LB, UB, relTol=1e-3):
    Phatss = np.atleast_2d(np.asarray(Phatss, float))
    M, d = Phatss.shape
    LB = np.asarray(LB, float)
    UB = np.asarray(UB, float)
    if LB.ndim == 1:
        LB = np.tile(LB, (M, 1))
    if UB.ndim == 1:
        UB = np.tile(UB, (M, 1))
    LB, UB = LB[:, :d].copy(), UB[:, :d].copy()
    flip = LB > UB
    LB[flip], UB[flip] = UB[flip].copy(), LB[flip].copy()
    hitLB = Phatss <= LB + relTol * np.abs(LB)
    hitUB = Phatss >= UB - relTol * np.abs(UB)
    valid = ~np.isnan(Phatss)
    hitLB &= valid
    hitUB &= valid
    nValid = valid.sum(0)
    fracLB = hitLB.sum(0) / np.maximum(nValid, 1)
    fracUB = hitUB.sum(0) / np.maximum(nValid, 1)
    pinned = np.all(LB == UB, axis=0)
    fracLB = fracLB.astype(float)
    fracUB = fracUB.astype(float)
    fracLB[pinned | (nValid == 0)] = np.nan
    fracUB[pinned | (nValid == 0)] = np.nan
    return fracLB, fracUB


def quantile_mcse(Phatss, B2=500, seed=20261030):
    """Monte Carlo SEs of the 2.5%/97.5% bootstrap quantiles and of the PII,
    by resampling the bootstrap estimates B2 times (own random stream)."""
    Phatss = np.atleast_2d(np.asarray(Phatss, float))
    M, d = Phatss.shape
    lo = np.full(d, np.nan)
    hi = np.full(d, np.nan)
    pii_ = np.full(d, np.nan)
    if M < 2 or B2 < 2:
        return lo, hi, pii_
    rs = np.random.RandomState(seed)
    idx = np.minimum(np.ceil(rs.random_sample((B2, M)).T * M).astype(int), M) - 1   # M x B2
    for k in range(d):
        x = Phatss[:, k]
        if np.all(np.isnan(x)) or np.nanmax(x) - np.nanmin(x) == 0:
            continue
        q = mquantile(x[idx], [0.025, 0.975], axis=0)         # 2 x B2
        lo[k] = np.std(q[0], ddof=1)
        hi[k] = np.std(q[1], ddof=1)
        ok = (q[0] > 0) & (q[1] > 0)
        if ok.sum() >= 2:
            pii_[k] = np.std(np.log10(q[1, ok]) - np.log10(q[0, ok]), ddof=1)
    return lo, hi, pii_
