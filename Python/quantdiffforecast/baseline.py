"""Flat and trend baseline forecasts (baselineForecast.m).

The baseline uses its own random stream, seeded by ``seed``, so it never
disturbs the random numbers of a run. MATLAB used RandStream('mt19937ar');
numpy's legacy ``RandomState`` is the same Mersenne Twister with the same
seeding, and integers are drawn as ceil(U * n) from its uniforms.
"""
from __future__ import annotations

import numpy as np


def _randi(rs, n, shape):
    u = rs.random_sample(size=(shape[1], shape[0])).T      # column-major fill, like MATLAB
    return np.minimum(np.ceil(u * n).astype(int), n) - 1   # 0-based


def baseline_forecast(data1, horizon, type="flat", Ndraws=6000, seed=1, k=None):
    """Return ((n + horizon) x Ndraws) draws and an info dict.

    'flat'  : last value plus sums of resampled, symmetrised past changes
    'trend' : log-linear fit to the last k (default min(7, n)) points,
              extrapolated, with resampled detrended log increments
    Rows 0..n-1 repeat the calibration data; draws are clipped at 0.
    """
    y = np.asarray(data1, dtype=float)[:, 1].ravel()
    n = y.size
    if n < 3:
        raise ValueError(f"baselineForecast:data - At least 3 calibration points are needed (got {n}).")
    if horizon < 1 or int(horizon) != horizon:
        raise ValueError("baselineForecast:horizon - horizon must be a positive integer.")
    if not np.all(np.isfinite(y)):
        raise ValueError("baselineForecast:data - Calibration data must be finite.")
    Ndraws = 2 * int(np.ceil(Ndraws / 2))
    half = Ndraws // 2
    rs = np.random.RandomState(seed)
    info = {"type": type, "lastValue": y[-1], "nDiffs": n - 1, "k": None, "slope": np.nan,
            "Ndraws": Ndraws, "seed": seed}
    t = type.lower()
    if t == "flat":
        dd = np.diff(y)
        sym = np.concatenate([dd, -dd])
        idx = _randi(rs, sym.size, (horizon, half))
        paths = np.cumsum(sym[idx], axis=0)
        fore = y[-1] + np.hstack([paths, -paths])
    elif t == "trend":
        if k is None:
            k = min(7, n)
        k = max(3, min(k, n))
        info["k"] = k
        tt = np.arange(1, n + 1, dtype=float)
        tk = tt[-k:]
        lk = np.log(y[-k:] + 1)
        X = np.column_stack([np.ones(k), tk])
        coef = np.linalg.lstsq(X, lk, rcond=None)[0]
        b = coef[1]
        info["slope"] = b
        mlog = coef[0] + b * (n + np.arange(1, horizon + 1))
        dl = np.diff(np.log(y + 1)) - b
        sym = np.concatenate([dl, -dl])
        idx = _randi(rs, sym.size, (horizon, half))
        paths = np.cumsum(sym[idx], axis=0)
        fore = np.exp(mlog[:, None] + np.hstack([paths, -paths])) - 1
    else:
        raise ValueError(f"baselineForecast:type - type must be 'flat' or 'trend' (got '{type}').")
    fore = np.maximum(fore, 0)
    draws = np.vstack([np.tile(y[:, None], (1, Ndraws)), fore])
    return draws, info
