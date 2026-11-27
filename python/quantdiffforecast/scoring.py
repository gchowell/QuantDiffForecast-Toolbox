"""Forecast scores: computeforecastperformance.m, computeWIS.m,
computeQuantiles.m, quantdiffForecastScores.m and computeScoresFromQuantiles.m.

MATLAB's ``quantile`` (Statistics Toolbox) places the sorted sample at
probabilities (i - 0.5)/n and interpolates linearly, clamping outside that
range. That is Hyndman-Fan type 5, numpy's ``method='hazen'``; every quantile
in this package goes through :func:`mquantile` so the numbers match MATLAB.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

WIS_ALPHAS = np.array([0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
ALPHA_QUANTILES = np.array([0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,
                            0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800,
                            0.850, 0.900, 0.950, 0.975, 0.990])
QUANTILE_NAMES = [f"Q_{q:.3f}" for q in ALPHA_QUANTILES]


def mquantile(x, p, axis=None):
    """MATLAB-compatible quantile (Hazen / type 5)."""
    return np.quantile(np.asarray(x, dtype=float), p, axis=axis, method="hazen")


def _wis_rows(D, y):
    """Vectorized WIS for rows of draws D (n x M) and observations y (n)."""
    K = WIS_ALPHAS.size
    lo = np.maximum(mquantile(D, WIS_ALPHAS / 2, axis=1).T, 0)        # n x K
    hi = np.maximum(mquantile(D, 1 - WIS_ALPHAS / 2, axis=1).T, 0)
    med = mquantile(D, 0.5, axis=1)
    yy = y[:, None]
    IS = (hi - lo) + (2 / WIS_ALPHAS) * (lo - yy) * (yy < lo) + (2 / WIS_ALPHAS) * (yy - hi) * (yy > hi)
    wis = (0.5 * np.abs(y - med) + np.sum((WIS_ALPHAS / 2) * IS, axis=1)) / (K + 0.5)
    return wis, med


def compute_forecast_performance(data1, datalatest, curves1, curves2, forecastingperiod):
    """computeforecastperformance: returns dict with calibration scalars
    (MAEc, MSEc, RMSEc, PIc, MISc) and per-horizon arrays [h value] for the
    forecast (cumulative over leads 1..h)."""
    tf2 = data1.shape[0]
    datac = datalatest[:tf2, 1]
    ym = np.median(curves1[:tf2, :], axis=1)
    out = {}
    out["RMSEC"] = np.sqrt(np.mean((datac - ym) ** 2))
    out["MSEC"] = np.mean((datac - ym) ** 2)
    out["MAEC"] = np.mean(np.abs(datac - ym))
    Lt = np.maximum(mquantile(curves2[:tf2, :], 0.025, axis=1), 0)
    Ut = np.maximum(mquantile(curves2[:tf2, :], 0.975, axis=1), 0)
    out["PIC"] = 100 * np.mean((datac >= Lt) & (datac <= Ut))
    out["MISC"] = np.mean((Ut - Lt) + (2 / 0.05) * (Lt - datac) * (datac < Lt)
                          + (2 / 0.05) * (datac - Ut) * (datac > Ut))
    for k in ("RMSEF", "MSEF", "MAEF", "PIF", "MISF"):
        out[k] = np.zeros((0, 2))
    if forecastingperiod > 0:
        if datalatest.shape[0] < tf2 + forecastingperiod:
            warnings.warn("forecasting period is too long to evaluate with available data")
            return out
        rows = {k: [] for k in ("RMSEF", "MSEF", "MAEF", "PIF", "MISF")}
        for h in range(1, forecastingperiod + 1):
            sl = slice(tf2, tf2 + h)
            yf = np.median(curves1[sl, :], axis=1)
            dataf = datalatest[sl, 1]
            L = np.maximum(mquantile(curves2[sl, :], 0.025, axis=1), 0)
            U = np.maximum(mquantile(curves2[sl, :], 0.975, axis=1), 0)
            rows["RMSEF"].append([h, np.sqrt(np.mean((dataf - yf) ** 2))])
            rows["MSEF"].append([h, np.mean((dataf - yf) ** 2)])
            rows["MAEF"].append([h, np.mean(np.abs(dataf - yf))])
            rows["PIF"].append([h, 100 * np.mean((dataf >= L) & (dataf <= U))])
            rows["MISF"].append([h, np.mean((U - L) + (2 / 0.05) * (L - dataf) * (dataf < L)
                                             + (2 / 0.05) * (dataf - U) * (dataf > U))])
        for k, v in rows.items():
            out[k] = np.array(v)
    return out


def compute_wis(data1, datalatest, curves2, forecastingperiod):
    """computeWIS: (mean calibration WIS, array [h, mean WIS over leads 1..h])."""
    tf2 = data1.shape[0]
    datac = datalatest[:tf2, 1]
    wisc, _ = _wis_rows(curves2[:tf2, :], datac)
    WISC = float(np.mean(wisc))
    WISFS = np.zeros((0, 2))
    if forecastingperiod > 0:
        if datalatest.shape[0] < tf2 + forecastingperiod:
            warnings.warn("forecasting period is too long to evaluate with available data")
            return WISC, WISFS
        wisf, _ = _wis_rows(curves2[tf2:tf2 + forecastingperiod, :],
                            datalatest[tf2:tf2 + forecastingperiod, 1])
        h = np.arange(1, forecastingperiod + 1)
        WISFS = np.column_stack([h, np.cumsum(wisf) / h])
    return WISC, WISFS


def compute_quantiles(data1, curves2, forecastingperiod):
    """computeQuantiles: 23 quantiles (clipped at 0) for calibration and forecast rows."""
    tf2 = data1.shape[0]
    qc = np.maximum(mquantile(curves2[:tf2, :], ALPHA_QUANTILES, axis=1).T, 0)
    qf = np.zeros((forecastingperiod, ALPHA_QUANTILES.size))
    if forecastingperiod > 0:
        qf = np.maximum(mquantile(curves2[tf2:tf2 + forecastingperiod, :], ALPHA_QUANTILES, axis=1).T, 0)
    return qc, qf, ALPHA_QUANTILES


# ---------------------------------------------------------------------------
# per-horizon scores (quantdiffForecastScores)
# ---------------------------------------------------------------------------
def _local_score(D, y):
    K = WIS_ALPHAS.size
    levels = np.concatenate([WIS_ALPHAS / 2, 1 - WIS_ALPHAS / 2, [0.5, 0.25, 0.75, 0.025, 0.975]])
    q = mquantile(D, levels, axis=1).T                     # n x levels
    L = np.maximum(q[:, :K], 0)
    U = np.maximum(q[:, K:2 * K], 0)
    med = q[:, 2 * K]
    yy = y[:, None]
    with np.errstate(invalid="ignore"):
        IS = (U - L) + (2 / WIS_ALPHAS) * (L - yy) * (yy < L) + (2 / WIS_ALPHAS) * (yy - U) * (yy > U)
        wis = (0.5 * np.abs(y - med) + np.sum((WIS_ALPHAS / 2) * IS, axis=1)) / (K + 0.5)
        c50 = ((y >= np.maximum(q[:, 2 * K + 1], 0)) & (y <= np.maximum(q[:, 2 * K + 2], 0))).astype(float)
        c95 = ((y >= np.maximum(q[:, 2 * K + 3], 0)) & (y <= np.maximum(q[:, 2 * K + 4], 0))).astype(float)
    wis[~np.isfinite(y)] = np.nan
    return wis, med, c50, c95


def _aggregate(F, k, nB):
    if k == 1:
        return F
    M = F.shape[1]
    return F.reshape(nB, k, M).sum(axis=1)


def forecast_scores(draws, datalatest, nCalib, horizon, aggregate=1, baseline_draws=None):
    """quantdiffForecastScores. Returns dict with DataFrames 'perHorizon' and
    'cumulative' and a dict 'summary'."""
    draws = np.atleast_2d(np.asarray(draws, dtype=float))
    datalatest = np.atleast_2d(np.asarray(datalatest, dtype=float))
    k = int(aggregate)
    if draws.shape[0] < nCalib + horizon:
        raise ValueError(f"quantdiffForecastScores:rows - draws has {draws.shape[0]} rows; "
                         f"{nCalib} calibration rows + {horizon} forecast steps are needed.")
    nB = horizon // k
    if nB < 1:
        raise ValueError(f"quantdiffForecastScores:aggregate - 'Aggregate' = {k} is longer "
                         f"than the horizon ({horizon} steps).")
    D = _aggregate(draws[nCalib:nCalib + nB * k, :], k, nB)
    nObs = datalatest.shape[0]
    y = np.full(nB, np.nan)
    tFrom = np.full(nB, np.nan)
    tTo = np.full(nB, np.nan)
    stepFrom = np.arange(nB) * k + 1
    stepTo = stepFrom + k - 1
    for b in range(nB):
        r0 = nCalib + stepFrom[b] - 1
        r1 = nCalib + stepTo[b] - 1
        if r0 < nObs:
            tFrom[b] = datalatest[r0, 0]
        if r1 < nObs:
            tTo[b] = datalatest[r1, 0]
            y[b] = datalatest[r0:r1 + 1, 1].sum()
    wis, med, c50, c95 = _local_score(D, y)
    with np.errstate(invalid="ignore"):
        wisLog, _, _, _ = _local_score(np.log1p(np.maximum(D, 0)), np.log1p(np.maximum(y, 0)))
    absErr = np.abs(y - med)
    sqErr = (y - med) ** 2
    scored = np.isfinite(y)
    c50[~scored] = np.nan
    c95[~scored] = np.nan
    wisLog[~scored] = np.nan
    per = pd.DataFrame({"horizon": np.arange(1, nB + 1), "step_from": stepFrom, "step_to": stepTo,
                        "time_from": tFrom, "time_to": tTo, "observed": y, "median": med,
                        "abs_error": absErr, "sq_error": sqErr, "cover50": c50, "cover95": c95,
                        "WIS": wis, "WIS_log": wisLog})
    if baseline_draws is not None:
        bd = np.asarray(baseline_draws, dtype=float)
        if bd.shape[0] < nCalib + horizon:
            raise ValueError("quantdiffForecastScores:baseline - BaselineDraws too short.")
        wisBase, _, _, _ = _local_score(_aggregate(bd[nCalib:nCalib + nB * k, :], k, nB), y)
        per["WIS_baseline"] = wisBase
        per["relWIS"] = wis / wisBase
    h = np.arange(1, nB + 1)
    cum = pd.DataFrame({"horizon": h, "MAE": np.cumsum(absErr) / h, "MSE": np.cumsum(sqErr) / h,
                        "coverage50": 100 * np.cumsum(c50) / h, "coverage95": 100 * np.cumsum(c95) / h,
                        "WIS": np.cumsum(wis) / h, "WIS_log": np.cumsum(wisLog) / h})

    def m(a):
        return float(np.mean(a[scored])) if scored.any() else np.nan
    summary = {"nScored": int(scored.sum()), "MAE": m(absErr), "MSE": m(sqErr),
               "RMSE": float(np.sqrt(m(sqErr))), "coverage50": 100 * m(c50),
               "coverage95": 100 * m(c95), "WIS": m(wis), "WIS_log": m(wisLog)}
    if baseline_draws is not None:
        both = scored & np.isfinite(per["WIS_baseline"].to_numpy())
        wb = float(np.mean(per["WIS_baseline"].to_numpy()[both])) if both.any() else np.nan
        summary["WIS_baseline"] = wb
        summary["relWIS"] = float(np.mean(wis[both]) / wb) if both.any() else np.nan
    return {"perHorizon": per, "cumulative": cum, "summary": summary}


def scores_from_quantiles(Q, levels, datalatest, nCalib, horizon):
    """computeScoresFromQuantiles: [MAE MSE coverage95 WIS] for calibration and
    the running means per lead (rows [h MAE MSE coverage WIS])."""
    levels = np.asarray(levels, dtype=float)

    def col(p):
        j = np.where(np.abs(levels - p) < 1e-9)[0]
        if j.size == 0:
            raise ValueError(f"The quantile level {p:g} is required but is not among the supplied levels.")
        return j[0]
    iMed = col(0.5)
    iLo = [col(a / 2) for a in WIS_ALPHAS]
    iHi = [col(1 - a / 2) for a in WIS_ALPHAS]
    iL95, iU95 = col(0.025), col(0.975)
    y = np.asarray(datalatest, dtype=float)[:, 1]
    nObs = y.size
    nRows = min(nObs, nCalib + horizon)
    Qr = Q[:nRows]
    yr = y[:nRows]
    m = Qr[:, iMed]
    absErr = np.abs(yr - m)
    sqErr = (yr - m) ** 2
    inside = ((yr >= np.maximum(Qr[:, iL95], 0)) & (yr <= np.maximum(Qr[:, iU95], 0))).astype(float)
    L = np.maximum(Qr[:, iLo], 0)
    U = np.maximum(Qr[:, iHi], 0)
    yy = yr[:, None]
    IS = (U - L) + (2 / WIS_ALPHAS) * (L - yy) * (yy < L) + (2 / WIS_ALPHAS) * (yy - U) * (yy > U)
    wis = (0.5 * absErr + np.sum((WIS_ALPHAS / 2) * IS, axis=1)) / (WIS_ALPHAS.size + 0.5)
    c = slice(0, nCalib)
    perfC = np.array([absErr[c].mean(), sqErr[c].mean(), 100 * inside[c].mean(), wis[c].mean()])
    perfF = None
    if horizon > 0 and nObs >= nCalib + horizon:
        perfF = np.zeros((horizon, 5))
        for h in range(1, horizon + 1):
            f = slice(nCalib, nCalib + h)
            perfF[h - 1] = [h, absErr[f].mean(), sqErr[f].mean(), 100 * inside[f].mean(), wis[f].mean()]
    return perfC, perfF
