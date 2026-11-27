"""Model comparison and forecast ensembles.

Python port of quantdiffCompareForecasts.m, ensembleWeights.m,
getensemble_linearpool_from_curves.m and run_ensemble.m.

Like the MATLAB functions, these read the saved forecast files
(Forecast-ODEModel-*.npz here) of each model, so models are compared on
exactly the predictive draws that were scored when they were run.
"""
from __future__ import annotations

import os
import re
import warnings

import numpy as np
import pandas as pd

from .baseline import baseline_forecast
from .config import Options
from .runners import forecast_npz_name, load_forecast, run_forecasting
from .scoring import forecast_scores, scores_from_quantiles

LEVELS_Q = np.concatenate([[0.01, 0.025], np.round(np.arange(0.05, 0.951, 0.05), 2), [0.975, 0.99]])


def _forecast_path(o: Options, fit_pos, tstart, tend, W, horizon):
    return os.path.join(o.outdir(), forecast_npz_name(
        o.cadfilename1, o.model.name, o.vars.fit_index[fit_pos], o.params.fixI0, o.method1,
        o.effective_dist1(), tstart, tend, W, horizon))


def compare_forecasts(options_list, windows, horizon, labels=None, tend=None, run=False, seed=1,
                      aggregate=1, baseline="flat", fit_index=0, save=True, name="comparison",
                      n_jobs=1, verbose=True):
    """quantdiffCompareForecasts.

    options_list  list of Options objects (one per model)
    windows       rows [first calibration row (1-based), number of calibration rows]
    fit_index     position within vars.fit_index of the series to score (0-based)
    Returns a dict with DataFrames perHorizon, byHorizon and summary.
    """
    if isinstance(options_list, Options):
        options_list = [options_list]
    K = len(options_list)
    windows = np.atleast_2d(np.asarray(windows))
    if windows.shape[1] != 2 or np.any(windows < 1) or np.any(windows != np.floor(windows)):
        raise ValueError("quantdiffCompareForecasts:windows - windows must have one row per "
                         "forecast origin: [tstart windowsize], positive integers.")
    windows = windows.astype(int)
    nW = windows.shape[0]
    labels = list(labels) if labels else [o.name or o.model.name for o in options_list]
    if len(labels) != K:
        raise ValueError("'labels' must have one entry per options object.")
    if tend is None:
        tends = windows[:, 0]
    else:
        tends = np.full(nW, tend) if np.isscalar(tend) else np.asarray(tend)
    baseline = baseline.lower()
    if baseline not in ("flat", "trend", "none"):
        raise ValueError("'baseline' must be 'flat', 'trend' or 'none'.")
    nH = horizon // aggregate
    if nH < 1:
        raise ValueError(f"'aggregate' = {aggregate} is longer than the horizon ({horizon} steps).")

    shape = (nW, nH, K)
    AE, SE, C50, C95, WIS, WISL, MED = (np.full(shape, np.nan) for _ in range(7))
    OBS = np.full((nW, nH), np.nan)
    WISB = np.full((nW, nH), np.nan)
    ORIGIN = np.full(nW, np.nan)
    files = [[None] * K for _ in range(nW)]
    for w in range(nW):
        tstart, W = int(windows[w, 0]), int(windows[w, 1])
        obsRef, baseDraws = None, None
        for k, o in enumerate(options_list):
            fpath = _forecast_path(o, fit_index, tstart, int(tends[w]), W, horizon)
            if not os.path.isfile(fpath) and run and tends[w] == tstart:
                if verbose:
                    print(f"compare_forecasts: running {labels[k]} for rows {tstart}-"
                          f"{tstart + W - 1}, horizon {horizon} ...")
                run_forecasting(o, tstart, tstart, W, horizon, seed=seed, n_jobs=n_jobs,
                                verbose=False, plot=False)
            if not os.path.isfile(fpath):
                raise FileNotFoundError(
                    f"quantdiffCompareForecasts:missing - No result file for {labels[k]}, "
                    f"calibration rows {tstart}-{tstart + W - 1}, horizon {horizon}:\n  {fpath}\n"
                    "Run run_forecasting first, or call this function with run=True.")
            F = load_forecast(fpath)
            files[w][k] = os.path.basename(fpath)
            nCalib = F["data1"].shape[0]
            obs = F["datalatest"][:, [0, fit_index + 1]]
            calib = F["data1"][:, [0, fit_index + 1]]
            nCmp = min(obs.shape[0], nCalib + horizon)
            if k == 0:
                obsRef = obs[:nCmp]
                ORIGIN[w] = obs[nCalib - 1, 0]
                if baseline != "none":
                    baseDraws, _ = baseline_forecast(calib, horizon, baseline,
                                                     F["forecast2"].shape[1], seed)
            elif nCalib != W or not np.array_equal(obs[:nCmp], obsRef):
                raise ValueError(f"quantdiffCompareForecasts:data - {labels[k]} was fitted to "
                                 f"different data than {labels[0]}.")
            S = forecast_scores(F["forecast2"], obs, nCalib, horizon, aggregate=aggregate,
                                baseline_draws=baseDraws)
            T = S["perHorizon"]
            if baseDraws is not None:
                WISB[w] = T["WIS_baseline"]
            OBS[w] = T["observed"]
            MED[w, :, k], AE[w, :, k], SE[w, :, k] = T["median"], T["abs_error"], T["sq_error"]
            C50[w, :, k], C95[w, :, k] = T["cover50"], T["cover95"]
            WIS[w, :, k], WISL[w, :, k] = T["WIS"], T["WIS_log"]

    rows = []
    for k in range(K):
        for w in range(nW):
            for h in range(nH):
                rows.append([labels[k], windows[w, 0], windows[w, 1], ORIGIN[w], h + 1, OBS[w, h],
                             MED[w, h, k], AE[w, h, k], SE[w, h, k], C50[w, h, k], C95[w, h, k],
                             WIS[w, h, k], WISL[w, h, k], WISB[w, h]])
    per = pd.DataFrame(rows, columns=["model", "tstart", "windowsize", "origin_time", "horizon",
                                      "observed", "median", "abs_error", "sq_error", "cover50",
                                      "cover95", "WIS", "WIS_log", "WIS_baseline"])
    per["relWIS"] = per["WIS"] / per["WIS_baseline"]

    common = np.all(np.isfinite(WIS), axis=2)
    if baseline != "none":
        common &= np.isfinite(WISB)

    def means(ae, se, c50, c95, wis, wisl, wisb, use):
        use = use.ravel().astype(bool)
        n = int(use.sum())
        if n == 0:
            return [np.nan] * 8 + [0]
        ae, se, c50, c95, wis, wisl, wisb = (a.ravel()[use] for a in (ae, se, c50, c95, wis, wisl, wisb))
        wb = np.mean(wisb)
        return [ae.mean(), np.sqrt(se.mean()), 100 * c50.mean(), 100 * c95.mean(), wis.mean(),
                wisl.mean(), wb, wis.mean() / wb, n]
    cols = ["MAE", "RMSE", "coverage50", "coverage95", "WIS", "WIS_log", "WIS_baseline", "relWIS"]
    by = [[labels[k], h + 1] + means(AE[:, h, k], SE[:, h, k], C50[:, h, k], C95[:, h, k],
                                      WIS[:, h, k], WISL[:, h, k], WISB[:, h], common[:, h])
          for k in range(K) for h in range(nH)]
    byH = pd.DataFrame(by, columns=["model", "horizon"] + cols + ["nWindows"])
    sm = [[labels[k]] + means(AE[:, :, k], SE[:, :, k], C50[:, :, k], C95[:, :, k], WIS[:, :, k],
                              WISL[:, :, k], WISB, common) for k in range(K)]
    summary = pd.DataFrame(sm, columns=["model"] + cols + ["nTargets"])
    R = {"perHorizon": per, "byHorizon": byH, "summary": summary, "windows": windows,
         "horizon": horizon, "aggregate": aggregate, "baseline": baseline, "files": files}
    if save:
        outdir = options_list[0].outdir()
        stem = os.path.join(outdir, "scores-comparison-" + re.sub(r'[\s/\\:*?"<>|]+', "_", name))
        per.to_csv(stem + "-perhorizon.csv", index=False)
        byH.to_csv(stem + "-byhorizon.csv", index=False)
        summary.to_csv(stem + "-summary.csv", index=False)
        R["savedTo"] = stem
    if verbose:
        print(f"\nForecast comparison: {nW} window(s), horizon {horizon}, baseline: {baseline}. "
              f"Means over the {int(common.sum())} target(s) scored for every model.")
        print(summary.to_string(index=False))
    return R


# ---------------------------------------------------------------------------
# ensembles
# ---------------------------------------------------------------------------
def ensemble_weights(weight_type1, AICcs=None, WISCs=None):
    """-1 equal; 0/1 Akaike; 2 inverse calibration WIS; 3 softmax calibration WIS."""
    A = np.asarray([] if AICcs is None else AICcs, float).ravel()
    Wc = np.asarray([] if WISCs is None else WISCs, float).ravel()
    if weight_type1 == -1:
        K = max(A.size, Wc.size)
        if K == 0:
            raise ValueError("ensembleWeights:K - Give AICcs or WISCs so the number of models is known.")
        vals, rule = np.ones(K), "equal"
    elif weight_type1 in (0, 1):
        if A.size == 0 or np.any(~np.isfinite(A)):
            raise ValueError("ensembleWeights:AICc - Akaike weights need finite AICc values.")
        vals, rule = np.exp(-(A - A.min()) / 2), "Akaike"
    elif weight_type1 in (2, 3):
        if Wc.size == 0:
            raise ValueError("ensembleWeights:WIS - WIS-based weights need the calibration WIS.")
        if np.any(~np.isfinite(Wc) | (Wc <= 0)):
            warnings.warn("ensembleWeights:wis - Non-finite or non-positive calibration WIS; "
                          "using equal weights.")
            vals, rule = np.ones(Wc.size), "equal (WIS fallback)"
        elif weight_type1 == 2:
            vals, rule = 1 / Wc, "inverse WIS"
        else:
            vals, rule = np.exp(-Wc / np.median(Wc)), "softmax WIS"
    else:
        raise ValueError(f"ensembleWeights:type - weight_type1 must be -1, 0, 1, 2 or 3 (got {weight_type1}).")
    return vals / vals.sum(), {"rule": rule, "raw": vals, "weight_type1": weight_type1}


def _weighted_quantiles(values, weights, probs):
    """Exact weighted empirical quantile inf{x : F(x) >= p}; values T x n."""
    order = np.argsort(values, axis=1, kind="stable")
    sv = np.take_along_axis(values, order, axis=1)
    cw = np.cumsum(weights[order] / weights.sum(), axis=1)
    cw[:, -1] = 1.0
    out = np.empty((values.shape[0], len(probs)))
    for j, p in enumerate(probs):
        if p <= 0:
            idx = np.zeros(values.shape[0], int)
        elif p >= 1:
            idx = np.full(values.shape[0], values.shape[1] - 1)
        else:
            idx = np.argmax(cw >= p, axis=1)
        out[:, j] = sv[np.arange(values.shape[0]), idx]
    return out


def linear_pool(weights, curves, probs_out=(0.5, 0.025, 0.975), clamp_nonneg=True, Nens=None, seed=1):
    """getensemble_linearpool_from_curves: weighted linear pool of the members'
    predictive draws.

    Returns (Qens, forecast1, curvesEns): quantiles at probs_out, the
    [median, 2.5%, 97.5%] summary, and Nens whole trajectories sampled from
    the members in proportion to the weights.
    """
    w = np.asarray(weights, float).ravel()
    K = len(curves)
    if w.size != K:
        raise ValueError("weights length must equal the number of model matrices.")
    if w.sum() <= 0 or abs(w.sum() - 1) > 1e-8:
        raise ValueError("weights must sum to one (tolerance 1e-8).")
    mats = [np.maximum(np.asarray(c, float), 0) if clamp_nonneg else np.asarray(c, float) for c in curves]
    T, M = mats[0].shape
    if any(m.shape != (T, M) for m in mats):
        raise ValueError("All model matrices must have the same size T x M.")
    active = np.where(w > 0)[0]
    pooled = np.hstack([mats[i] for i in active])
    pw = np.concatenate([np.full(M, w[i] / M) for i in active])
    probs_out = np.asarray(probs_out, float)
    Qens = _weighted_quantiles(pooled, pw, probs_out)
    forecast1 = _weighted_quantiles(pooled, pw, [0.5, 0.025, 0.975])
    Nens = M if Nens is None else int(Nens)
    curvesEns = None
    if Nens > 0:
        rs = np.random.RandomState(int(seed))
        aw = w[active] / w[active].sum()
        cum = np.cumsum(aw)
        cum[-1] = 1.0
        u = rs.random_sample(Nens)
        lab = np.searchsorted(cum, u, side="left") if active.size > 1 else np.zeros(Nens, int)
        curvesEns = np.zeros((T, Nens))
        for a, i in enumerate(active):
            loc = np.where(lab == a)[0]
            if loc.size:
                chosen = np.minimum(np.ceil(rs.random_sample(loc.size) * M).astype(int), M) - 1
                curvesEns[:, loc] = mats[i][:, chosen]
    return Qens, forecast1, curvesEns


def run_ensemble(options_list, model_names, weight_type1, calibrationperiods, forecastingperiod,
                 name_ensemble, tstart=1, tend=1, Nens=None, seed=1, plot=True, fit_index=0,
                 verbose=True):
    """run_ensemble: pool the saved forecasts of several models.

    Each member must already have been run (run_forecasting) for the same
    window and horizon. Returns a dict with weights, member and ensemble scores
    and the relative WIS against the flat and trend baselines.
    """
    K = len(options_list)
    if len(model_names) != K:
        raise ValueError("model_names must have one entry per options object.")
    cps = np.atleast_1d(calibrationperiods).astype(int)
    nW = cps.size
    outdir = options_list[0].outdir()
    R = {"calibrationperiods": cps, "members": list(model_names), "weight_type1": weight_type1,
         "rule": "", "weights": np.full((nW, K), np.nan), "AICcs": np.full((nW, K), np.nan),
         "performanceC_models": np.full((nW, 4, K), np.nan),
         "performanceF_models": np.full((nW, 4, K), np.nan),
         "performanceC_ens": np.full((nW, 4), np.nan), "performanceF_ens": np.full((nW, 4), np.nan),
         "performanceF_ens_byhorizon": [None] * nW, "WIS_baselines": np.full((nW, 2), np.nan),
         "relWIS_ens": np.full((nW, 2), np.nan), "relWIS_models": np.full((nW, K, 2), np.nan),
         "forecast": [None] * nW}
    for wi, W in enumerate(cps):
        S, curves, AICc, WISC = [], [], np.full(K, np.nan), np.full(K, np.nan)
        for k, o in enumerate(options_list):
            fpath = _forecast_path(o, fit_index, tstart, tend, int(W), forecastingperiod)
            if not os.path.isfile(fpath):
                raise FileNotFoundError(f"run_ensemble:missing - Result file for model {k + 1} "
                                        f"({model_names[k]}), window {W} not found:\n  {fpath}\n"
                                        "Run run_forecasting for it first.")
            F = load_forecast(fpath)
            S.append(F)
            curves.append(np.maximum(F["forecast2"], 0))
            AICc[k] = F["AICcs"][-1, 1]
            WISC[k] = F["WISCSS"][-1]
        R["AICcs"][wi] = AICc
        Ms = [c.shape[1] for c in curves]
        if len({c.shape[0] for c in curves}) != 1:
            raise ValueError("run_ensemble:grid - Members have different forecast lengths.")
        if len(set(Ms)) != 1:
            warnings.warn(f"run_ensemble:draws - Members have different numbers of draws; "
                          f"truncating to {min(Ms)}.")
            curves = [c[:, :min(Ms)] for c in curves]
        weights1, winfo = ensemble_weights(weight_type1, AICc, WISC)
        R["weights"][wi] = weights1
        R["rule"] = winfo["rule"]
        NensW = curves[0].shape[1] if Nens is None else Nens
        Qens, forecast1, curvesEns = linear_pool(weights1, curves, LEVELS_Q, Nens=NensW, seed=seed)
        data1 = S[0]["data1"][:, [0, fit_index + 1]]
        datalatest = S[0]["datalatest"][:, [0, fit_index + 1]]
        for k in range(1, K):
            if not np.array_equal(S[k]["datalatest"][:, [0, fit_index + 1]], datalatest):
                raise ValueError(f"run_ensemble:data - Model {k + 1} ({model_names[k]}) was fitted "
                                 "to different data than model 1.")
        nCalib = data1.shape[0]
        perfC_e, perfF_e = scores_from_quantiles(Qens, LEVELS_Q, datalatest, nCalib, forecastingperiod)
        R["performanceC_ens"][wi] = perfC_e
        hasF = perfF_e is not None
        if hasF:
            R["performanceF_ens"][wi] = perfF_e[-1, 1:5]
        for k in range(K):
            Qk, _, _ = linear_pool([1.0], [curves[k]], LEVELS_Q, Nens=0)
            pc, pf = scores_from_quantiles(Qk, LEVELS_Q, datalatest, nCalib, forecastingperiod)
            R["performanceC_models"][wi, :, k] = pc
            if pf is not None:
                R["performanceF_models"][wi, :, k] = pf[-1, 1:5]
        o0 = options_list[0]
        sfx = (f"-model_name-{name_ensemble}-weight_type1-{weight_type1}-method-{o0.method1}-dist-"
               f"{o0.effective_dist1()}-tstart-{tstart}-tend-{tend}-calibrationperiod-{W}-horizon-"
               f"{forecastingperiod}-{o0.caddisease}-{o0.datatype}.csv")
        tv, dA = S[0]["timevect_all"], S[0]["data_all"][:, fit_index]
        nF = nCalib + forecastingperiod
        if tv.size < nF:
            pad = np.full(nF - tv.size, np.nan)
            tv, dA = np.concatenate([tv, pad]), np.concatenate([dA, pad])
        ftab = pd.DataFrame({"time": tv[:nF], "data": dA[:nF], "median": forecast1[:, 0],
                             "LB": forecast1[:, 1], "UB": forecast1[:, 2]})
        ftab.to_csv(os.path.join(outdir, "Forecast" + sfx), index=False)
        R["forecast"][wi] = ftab
        ctab = pd.DataFrame([[W, *perfC_e, *weights1]],
                            columns=["Calibrationperiod", "MAE", "MSE", "Coverage 95%PI", "WIS"]
                            + [f"w_{re.sub(r'[^0-9a-zA-Z_]', '_', n)}" for n in model_names])
        ctab.to_csv(os.path.join(outdir, "performance-calibration" + sfx), index=False)
        if hasF:
            calib = datalatest[:nCalib]
            WISb = {}
            for bt in ("flat", "trend"):
                bd, _ = baseline_forecast(calib, forecastingperiod, bt, NensW, seed)
                Qb, _, _ = linear_pool([1.0], [bd], LEVELS_Q, Nens=0)
                _, pfb = scores_from_quantiles(Qb, LEVELS_Q, datalatest, nCalib, forecastingperiod)
                WISb[bt] = pfb[:, 4]
            R["WIS_baselines"][wi] = [WISb["flat"][-1], WISb["trend"][-1]]
            R["relWIS_ens"][wi] = perfF_e[-1, 4] / R["WIS_baselines"][wi]
            for k in range(K):
                R["relWIS_models"][wi, k] = R["performanceF_models"][wi, 3, k] / R["WIS_baselines"][wi]
            perfF = np.column_stack([perfF_e, WISb["flat"], WISb["trend"],
                                     perfF_e[:, 4] / WISb["flat"], perfF_e[:, 4] / WISb["trend"]])
            ptab = pd.DataFrame(perfF, columns=["Horizon", "MAE", "MSE", "Coverage 95%PI", "WIS",
                                                "WIS_baseline_flat", "WIS_baseline_trend",
                                                "relWIS_flat", "relWIS_trend"])
            ptab.to_csv(os.path.join(outdir, "performance-forecasting" + sfx), index=False)
            R["performanceF_ens_byhorizon"][wi] = ptab
        if plot:
            from .plotting import plot_ensemble
            plot_ensemble(curves, list(model_names), weights1, forecast1,
                          f"{name_ensemble} ({winfo['rule']})", S[0]["timevect1"], S[0]["timevect2"],
                          datalatest, forecastingperiod,
                          path=os.path.join(outdir, f"fig-ensemble-{re.sub(r'[^0-9a-zA-Z_-]', '_', name_ensemble)}"
                                                    f"-calibrationperiod-{W}.png"))
        if verbose:
            print(f"Window {W}: weights {np.round(weights1, 3)} ({winfo['rule']}); ensemble "
                  f"calibration WIS {perfC_e[3]:.3g} vs members {np.round(WISC, 3)}")
            if hasF:
                print(f"           horizon {forecastingperiod} relative WIS vs flat / trend baseline: "
                      f"ensemble {R['relWIS_ens'][wi, 0]:.2f} / {R['relWIS_ens'][wi, 1]:.2f}")
    return R
