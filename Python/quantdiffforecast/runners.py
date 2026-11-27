"""Fitting and forecasting runners.

Python port of Run_Fit_ODEModel.m and Run_Forecasting_ODEModel.m. Both MATLAB
runners execute the same rolling-window algorithm (Run_Fit_ODEModel is the
forecasting runner with a zero horizon), so here they share one engine,
:func:`_run_windows`:

  for each window start i = tstart1 .. tend1
    1. fit the observed calibration rows (MultiStart, numstartpoints)
    2. AICc, and the noise model of the bootstrap (dist1 / factor1)
    3. B parametric-bootstrap datasets -> refit each one (1 start point)
       -> forecast curve per refit + 20 noisy draws per refit
    4. scores (MAE, MSE, 95% PI coverage, WIS, per-lead scores, baseline)
    5. parameter intervals, PII / bound-hit / Monte Carlo diagnostics
  then write the same CSV files as MATLAB (plus .npz instead of .mat).

What changed relative to MATLAB, and why
  * no globals: the options object carries method1, and the objective uses an
    explicit context, so the bootstrap refits can run in parallel (n_jobs);
  * random numbers come from one numpy Generator (seed=...); bootstrap
    datasets and noise draws are generated in the parent process in a fixed
    order, and each refit receives its own child seed, so results do not depend
    on n_jobs;
  * .mat files become .npz files with the same variable names;
  * figures are drawn with matplotlib (see plotting.py) and saved as PNG.
"""
from __future__ import annotations

import json
import os
import re
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from .baseline import baseline_forecast
from .config import (Options, QuantDiffForecastError, check_initial_states, check_negative_data,
                     load_data, make_window, validate_options)
from .core import aicc_count_adjustment, effective_n, get_aicc, get_nparams
from .identifiability import bound_hit_fractions, pii, quantile_mcse
from .noise import (add_error_structure, estimate_laplace_scale, estimate_pearson_dispersion,
                    negative_noise_share)
from .optimizer import fit_model
from .scoring import (QUANTILE_NAMES, compute_forecast_performance, compute_quantiles,
                      compute_wis, forecast_scores, mquantile)

DT = 1.0


# ---------------------------------------------------------------------------
# file names (identical to the MATLAB ones, .mat -> .npz)
# ---------------------------------------------------------------------------
def _num(x):
    return str(int(x)) if float(x).is_integer() else f"{x:g}"


def forecast_npz_name(cadfilename1, modelname, fit_index0, fixI0, method1, dist1, tstart, tend,
                      windowsize, horizon):
    """qdf_forecast_matname, with fit_index given 0-based (written 1-based)."""
    cad = cadfilename1 if cadfilename1.lower().endswith(".txt") else cadfilename1 + ".txt"
    return (f"Forecast-ODEModel-{cad}-model_name-{modelname}-fit_index-{int(fit_index0) + 1}"
            f"-fixI0-{_num(fixI0)}-method-{_num(method1)}-dist-{_num(dist1)}-tstart-{_num(tstart)}"
            f"-tend-{_num(tend)}-calibrationperiod-{_num(windowsize)}-forecastingperiod-{_num(horizon)}.npz")


def _safe_label(s):
    return re.sub(r'[\\/:*?"<>|]+', "", str(s))


# ---------------------------------------------------------------------------
# results containers
# ---------------------------------------------------------------------------
@dataclass
class SeriesResult:
    fit_index: int                     # 0-based state position
    data1: np.ndarray                  # [time, observed] calibration rows
    datalatest: np.ndarray             # [time, observed] from the window start onwards
    fit1: np.ndarray                   # bootstrap fitted curves (calibration rows)
    forecast1: np.ndarray              # bootstrap model trajectories (no noise)
    forecast2: np.ndarray              # predictive draws (with observation noise)
    timevect1: np.ndarray
    timevect2: np.ndarray
    forecast_table: pd.DataFrame       # time, data, median, LB, UB, model_median
    performance_calibration: pd.DataFrame
    performance_forecasting: Optional[pd.DataFrame] = None
    scores_by_horizon: Optional[pd.DataFrame] = None
    quantiles: Optional[pd.DataFrame] = None


@dataclass
class WindowResult:
    tstart: int
    P: np.ndarray                      # point estimate [params, X0, alpha, d]
    fval: float
    Phatss: np.ndarray                 # B x nz bootstrap estimates
    fvals: np.ndarray
    param_estims: np.ndarray           # nz x 3: estimate, 2.5%, 97.5%
    param_bootmedian: np.ndarray
    composite: np.ndarray              # [value, 2.5%, 97.5%]
    AICc: tuple
    noise: dict
    fit_diagnostics: dict
    series: list = field(default_factory=list)
    Ys: Optional[list] = None          # per state: (len(timevect2) x B) bootstrap states
    pii: dict = field(default_factory=dict)


@dataclass
class RunResult:
    options: Options
    tstart1: int
    tend1: int
    windowsize1: int
    forecastingperiod: int
    windows: list
    AICcs: pd.DataFrame
    parameters: pd.DataFrame
    composite: pd.DataFrame
    output_dir: str
    accepted_bootstrap: int = 0
    requested_bootstrap: int = 0

    # MATLAB-style accessors -------------------------------------------------
    @property
    def forecast_model1(self):
        return np.vstack([s.forecast1 for s in self.windows[-1].series])

    @property
    def forecast_model12(self):
        return np.vstack([s.forecast2 for s in self.windows[-1].series])

    @property
    def performanceC(self):
        return pd.concat([s.performance_calibration for w in self.windows for s in w.series],
                         ignore_index=True)

    @property
    def performanceF(self):
        tabs = [s.performance_forecasting for w in self.windows for s in w.series
                if s.performance_forecasting is not None]
        return pd.concat(tabs, ignore_index=True) if tabs else None


# ---------------------------------------------------------------------------
# bootstrap worker (module level so it can be sent to other processes)
# ---------------------------------------------------------------------------
def _boot_fit(data1_boot, P0, model, params, vars, horizon, method1, seed):
    rng = np.random.default_rng(seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = fit_model(data1_boot, P0, 1, DT, model, params, vars, horizon, method1, rng=rng)
    dg = r.diagnostics
    keep = {k: dg[k] for k in ("optimizerExitflag", "finalObjective", "converged", "LB", "UB",
                                "startAgreement", "scaledGradientInfNorm")}
    return r.P, r.fitcurve, r.forecastcurve, r.F2, r.fval, keep, r.timevect2


def _parallel_map(func, tasks, n_jobs):
    if n_jobs is None or n_jobs == 1 or len(tasks) < 2:
        return [func(*t) for t in tasks]
    try:
        from joblib import Parallel, delayed
    except ImportError:  # pragma: no cover
        warnings.warn("joblib is not installed: running the bootstrap serially.")
        return [func(*t) for t in tasks]
    return Parallel(n_jobs=n_jobs)(delayed(func)(*t) for t in tasks)


# ---------------------------------------------------------------------------
# the engine
# ---------------------------------------------------------------------------
def _run_windows(opts: Options, tstart1, tend1, windowsize1, forecastingperiod, getperformance,
                 caller, seed=None, rng=None, n_jobs=1, verbose=True, plot=None, save=True):
    validate_options(opts, caller)
    params, vars, model = opts.params, opts.vars, opts.model
    method1 = opts.method1
    dist1 = opts.effective_dist1()
    d = 1.0
    M = int(opts.B)
    nfit = len(vars.fit_index)
    caddisease, datatype = opts.caddisease, opts.datatype
    cad = opts.cadfilename_txt()
    outdir = opts.outdir()
    printscreen1 = opts.printscreen1 if plot is None else plot
    if rng is None:
        rng = np.random.default_rng(seed)

    data = load_data(opts)
    nrows = data.shape[0]
    numparams = get_nparams(method1, dist1, int(np.sum(params.fixed == 0)),
                            nfit * (params.fixI0 == 0))
    nEffWindow, nEffInfo = effective_n(windowsize1 * nfit, nfit, params.fixI0)
    if numparams >= nEffWindow - 1:
        raise QuantDiffForecastError(
            f"The calibration window has {nEffWindow} informative data points but {numparams} "
            f"quantities are estimated. AICc needs at least {numparams + 2}: increase the "
            "calibration period or fix more parameters.")
    if tend1 < tstart1:
        raise QuantDiffForecastError("tstart1 parameter should not exceed tend1")
    need = windowsize1 - 1 + forecastingperiod * (1 if getperformance else 0)
    if tend1 + need > nrows:
        tend_new = nrows - need
        if tend_new < tstart1:
            raise QuantDiffForecastError(
                f"No rolling window fits: rows through {tstart1 + need} are needed but the data "
                f"has {nrows} rows.")
        warnings.warn(f"{caller}:tend1 - tend1 = {tend1} exceeds the data ({nrows} rows); "
                      f"using tend1 = {tend_new}.")
        tend1 = tend_new
    check_negative_data(method1, data, range(tstart1, tend1 + windowsize1), caller)
    check_initial_states(vars, tend1, caller)

    nz = params.num + nfit + 2
    nslot = params.num + nfit + (1 if method1 in (3, 4) else 2 if method1 == 5 else 0)
    piiUser = np.full(nslot, np.nan)
    if params.PIIscale is not None:
        su = np.asarray(params.PIIscale, float).ravel()[:nslot]
        piiUser[:su.size] = su

    windows, AICcs, noiseLog = [], [], []
    accepted = 0
    for i in range(tstart1, tend1 + 1):
        i0 = i - 1
        rows = slice(i0, i0 + windowsize1)
        datac = data[rows, 1:]
        timevect1 = DT * data[rows, 0]
        data_all = data[i0:, 1:]
        timevect_all = DT * data[i0:, 0]
        if datac[0, 0] == 0:
            warnings.warn(f"window {i}: first data point in the time series is zero")
        if getperformance and data_all.shape[0] < windowsize1 + forecastingperiod:
            raise QuantDiffForecastError("Length of time series data is too short to evaluate "
                                         "the forecasting period indicated in forecastingperiod.")
        ydata = datac.flatten(order="F")
        I0 = datac[0, :]
        vars.window = make_window(data, i, windowsize1, DT)
        data1 = np.column_stack([data[rows, 0], datac])
        if params.fixI0:
            params0 = np.concatenate([params.initial, I0, [1, 1]])
        else:
            params0 = np.concatenate([params.initial, vars.initial[vars.fit_index], [1, 1]])

        if verbose:
            print(f"[{caller}] window {i}: fitting rows {i}-{i + windowsize1 - 1} "
                  f"(method {method1}, {opts.numstartpoints} start points)")
        fit = fit_model(data1, params0, opts.numstartpoints, DT, model, params, vars, 0, method1,
                        rng=rng)
        P_hat, fval, fitcurve_hat = fit.P, fit.fval, fit.fitcurve
        LBbox, UBbox = fit.diagnostics["LB"], fit.diagnostics["UB"]

        nEff, _ = effective_n(ydata.size, nfit, params.fixI0)
        adj = aicc_count_adjustment(method1, ydata, nfit, params.fixI0,
                                    P_hat[params.num + nfit], P_hat[params.num + nfit + 1])
        AICc, part1, part2, kpar = get_aicc(method1, dist1, int(np.sum(params.fixed == 0)),
                                            nfit * (params.fixI0 == 0), fval, nEff, adj)
        AICcs.append([i, AICc, part1, part2, kpar])

        # ---- noise model of the bootstrap -------------------------------
        dist1_boot = None
        noiseRaw = np.nan
        factor1 = 1.0
        if method1 == 0 and dist1 == 1:
            factor1 = 1.0
        elif method1 == 0 and dist1 == 2:
            factor1, dinfo = estimate_pearson_dispersion(ydata, fitcurve_hat, kpar - 1, nfit,
                                                         params.fixI0)
            noiseRaw = factor1
            if not np.isfinite(factor1):
                warnings.warn(f"QuantDiffForecast:dispersion - Window {i}: the dispersion cannot "
                              "be estimated; using Poisson bootstrap noise for this window.")
            if not np.isfinite(factor1) or factor1 <= 1:
                dist1_boot, factor1 = 1, 1.0
        elif method1 == 0 and dist1 == 0:
            var1 = np.sum((fitcurve_hat - ydata) ** 2) / (nEff - (kpar - 1))
            factor1 = float(np.sqrt(var1))
        elif method1 == 1:
            factor1 = 1.0
        elif method1 in (3, 4):
            factor1 = P_hat[params.num + nfit]
        elif method1 == 5:
            factor1 = P_hat[params.num + nfit]
            d = P_hat[params.num + nfit + 1]
        elif method1 == 6:
            factor1, linfo = estimate_laplace_scale(ydata, fitcurve_hat, kpar - 1, nfit, params.fixI0)
            noiseRaw = factor1
            if not np.isfinite(factor1):
                factor1 = np.sum(np.abs(fitcurve_hat - ydata)) / ydata.size
                warnings.warn(f"QuantDiffForecast:laplaceScale - Window {i}: using the "
                              f"uncorrected Laplace scale SAD/n = {factor1:g}.")
            if not factor1 > 0:
                raise QuantDiffForecastError(f"QuantDiffForecast:laplaceScale - Window {i}: the "
                                             "Laplace scale is 0 (an exact fit).")
        if dist1_boot is None:
            dist1_boot = dist1
        negShare, negInfo = negative_noise_share(fitcurve_hat, factor1, dist1_boot, nfit)
        if np.isfinite(negShare) and negShare > 0.05:
            warnings.warn(f"QuantDiffForecast:negativeNoise - Window {i}: about "
                          f"{100 * negShare:.0f}% of the bootstrap draws will be negative. "
                          "Consider dist1 = 2 with method1 = 0, or a count likelihood.")
        noise = {"tstart": i, "dist1": dist1, "dist1_boot": dist1_boot, "factor1": float(factor1),
                 "d": float(d), "rawEstimate": float(noiseRaw), "negativeShare": float(negShare)}
        noiseLog.append([i, dist1, dist1_boot, factor1, d, noiseRaw, negShare])

        # ---- parametric bootstrap -----------------------------------------
        if verbose:
            print(f"[{caller}] window {i}: bootstrapping {M} datasets"
                  + (f" on {n_jobs} workers" if n_jobs not in (None, 1) else ""))
        cum = np.cumsum(fitcurve_hat)
        tasks = []
        for j in range(M):
            sim = add_error_structure(cum, 1, dist1_boot, factor1, d, rng)[:, 0]
            data1_boot = np.column_stack([data[rows, 0],
                                          sim.reshape(windowsize1, nfit, order="F")])
            data1_boot[0, 1:] = I0
            tasks.append((data1_boot, P_hat.copy(), model, params, vars, forecastingperiod,
                          method1, int(rng.integers(2 ** 63 - 1))))
        outs = _parallel_map(_boot_fit, tasks, n_jobs)

        Phatss = np.zeros((M, nz))
        fvals = np.zeros(M)
        LBboot = np.full((M, nz), np.nan)
        UBboot = np.full((M, nz), np.nan)
        fit_model1, forecast_model1, forecast_model12, F2s, boot_diag = [], [], [], [], []
        timevect2 = None
        for j, (Pj, fc, fcast, F2, fv, dg, tv2) in enumerate(outs):
            Phatss[j] = Pj
            fvals[j] = fv
            LBboot[j], UBboot[j] = dg["LB"][:nz], dg["UB"][:nz]
            fit_model1.append(fc)
            forecast_model1.append(fcast)
            F2s.append(F2)
            boot_diag.append(dg)
            timevect2 = tv2
            if (method1 == 0 and dist1_boot in (0, 2)) or (method1 == 6 and dist1_boot == 6):
                forecast_model12.append(add_error_structure(np.cumsum(fcast), 20, dist1_boot,
                                                            factor1, 0, rng))
            else:
                forecast_model12.append(add_error_structure(np.cumsum(fcast), 20, dist1_boot,
                                                            Pj[-2], Pj[-1], rng))
            if dg["optimizerExitflag"] > 0 and np.isfinite(dg["finalObjective"]):
                accepted += 1
        fit_model1 = np.column_stack(fit_model1)
        forecast_model1 = np.column_stack(forecast_model1)
        forecast_model12 = np.hstack(forecast_model12)
        Ys = [np.column_stack([F[:, s] for F in F2s]) for s in range(vars.num)]

        data1 = np.column_stack([timevect1, datac])          # observed, time scaled by DT
        datalatest = data[i0:, :]

        if save:
            np.savez_compressed(
                os.path.join(outdir, f"bootstraps-ODEModel-{cad}-model_name-{model.name}-fixI0-"
                             f"{params.fixI0}-method-{method1}-dist-{dist1}-tstart-{i}-tend-{tend1}"
                             f"-calibrationperiod-{windowsize1}-forecastingperiod-{forecastingperiod}.npz"),
                Phatss_model1=Phatss, fvals_model1=fvals, LBbox=LBbox, UBbox=UBbox,
                LBboot=LBboot, UBboot=UBboot, noiseLog=np.array(noiseLog, dtype=float))

        # ---- scores per fitted series -------------------------------------
        series = []
        nC = windowsize1
        nF = windowsize1 + forecastingperiod
        for j in range(nfit):
            fit1 = fit_model1[j * nC:(j + 1) * nC, :]
            forecast1 = forecast_model1[j * nF:(j + 1) * nF, :]
            forecast2 = forecast_model12[j * nF:(j + 1) * nF, :]
            d1 = data1[:, [0, j + 1]]
            dl = datalatest[:, [0, j + 1]]
            perf = compute_forecast_performance(d1, dl, forecast2, forecast2, forecastingperiod)
            WISC, WISFS = compute_wis(d1, dl, forecast2, forecastingperiod)
            sfx = (f"model_name-{model.name}-fit_index-{vars.fit_index[j] + 1}-fixI0-{params.fixI0}"
                   f"-method-{method1}-dist-{dist1}-tstart-{i}-tend-{tend1}-calibrationperiod-"
                   f"{windowsize1}-horizon-{forecastingperiod}-{caddisease}-{datatype}.csv")
            perfF_tab, scoresH = None, None
            if forecastingperiod > 0 and perf["MAEF"].shape[0] > 0 and getperformance:
                perfF_tab = pd.DataFrame({
                    "forecasting_horizon": np.arange(1, forecastingperiod + 1),
                    "MAE": perf["MAEF"][:, 1], "MSE": perf["MSEF"][:, 1],
                    "Coverage 95%PI": perf["PIF"][:, 1], "WIS": WISFS[:, 1]})
                try:
                    bd, _ = baseline_forecast(d1, forecastingperiod, "flat", forecast2.shape[1], 1)
                except Exception:
                    bd = None
                scoresH = forecast_scores(forecast2, dl, nC, forecastingperiod,
                                          baseline_draws=bd)["perHorizon"]
                if save:
                    perfF_tab.to_csv(os.path.join(outdir, "performance-forecasting-" + sfx), index=False)
                    scoresH.to_csv(os.path.join(outdir, "scores-byhorizon-" + sfx), index=False)
            perfC_tab = pd.DataFrame([[i, windowsize1, perf["MAEC"], perf["MSEC"], perf["PIC"], WISC]],
                                     columns=["time", "calibration_period", "MAE", "MSE",
                                              "Coverage 95%PI", "WIS"])
            UB1 = np.maximum(mquantile(forecast2, 0.975, axis=1), 0)
            LB1 = np.maximum(mquantile(forecast2, 0.025, axis=1), 0)
            med1 = np.maximum(np.median(forecast2, axis=1), 0)
            mmed1 = np.maximum(np.median(forecast1, axis=1), 0)
            if data_all.shape[0] >= nF:
                tcol, dcol = timevect_all[:nF], data_all[:nF, j]
            else:
                pad = np.full(nF - data_all.shape[0], np.nan)
                tcol = np.concatenate([timevect_all, pad])
                dcol = np.concatenate([data_all[:, j], pad])
            ftab = pd.DataFrame({"time": tcol, "data": dcol, "median": med1, "LB": LB1, "UB": UB1,
                                 "model_median": mmed1})
            qc, qf, _ = compute_quantiles(d1, forecast2, forecastingperiod)
            qtab = pd.DataFrame(np.vstack([qc, qf]), columns=QUANTILE_NAMES)
            if save:
                perfC_tab.to_csv(os.path.join(outdir, "performance-calibration-" + sfx), index=False)
                ftab.to_csv(os.path.join(outdir, "Forecast-" + sfx), index=False)
                qtab.to_csv(os.path.join(outdir, "quantile-" + sfx), index=False)
                np.savez_compressed(
                    os.path.join(outdir, forecast_npz_name(cad, model.name, vars.fit_index[j],
                                                           params.fixI0, method1, dist1, i, tend1,
                                                           windowsize1, forecastingperiod)),
                    forecast1=forecast1, forecast2=forecast2, fit1=fit1, data1=data1,
                    datalatest=datalatest, timevect1=timevect1, timevect2=timevect2,
                    timevect_all=timevect_all, data_all=data_all,
                    forecastingperiod=forecastingperiod, AICcs=np.array(AICcs, float),
                    MAECSS=[perf["MAEC"]], MSECSS=[perf["MSEC"]], PICSS=[perf["PIC"]],
                    WISCSS=[WISC],
                    MAEFSS=perf["MAEF"][-1:, 1] if perf["MAEF"].size else np.zeros(0),
                    MSEFSS=perf["MSEF"][-1:, 1] if perf["MSEF"].size else np.zeros(0),
                    PIFSS=perf["PIF"][-1:, 1] if perf["PIF"].size else np.zeros(0),
                    WISFSS=WISFS[-1:, 1] if WISFS.size else np.zeros(0),
                    P_model1d=P_hat, Phatss_model1=Phatss,
                    meta=json.dumps({"model": model.name, "options": opts.name, "method1": method1,
                                     "dist1": dist1, "fit_index0": int(vars.fit_index[j]),
                                     "caddisease": caddisease, "datatype": datatype}))
            series.append(SeriesResult(
                fit_index=int(vars.fit_index[j]), data1=d1, datalatest=dl, fit1=fit1,
                forecast1=forecast1, forecast2=forecast2, timevect1=timevect1, timevect2=timevect2,
                forecast_table=ftab, performance_calibration=perfC_tab,
                performance_forecasting=perfF_tab, scores_by_horizon=scoresH, quantiles=qtab))

        # ---- parameter summaries -------------------------------------------
        q025 = mquantile(Phatss, 0.025, axis=0)
        q975 = mquantile(Phatss, 0.975, axis=0)
        param_estims = np.column_stack([P_hat, q025, q975])
        bootmed = np.median(Phatss, axis=0)
        outside = (P_hat < q025) | (P_hat > q975)
        if outside.any():
            zlab = params.label + [f"X0({k + 1})" for k in vars.fit_index] + ["alpha", "d"]
            warnings.warn(f"QuantDiffForecast:estimateOutsideCI - Window {i}: the fitted value of "
                          f"{', '.join(np.array(zlab)[outside])} lies outside its bootstrap 95% "
                          "interval (bootstrap bias; often an estimate on a bound or a ridge).")
        PIIv, rangeFrac, scale, PIIbox, ratio = pii(q025[:nslot], q975[:nslot], LBbox[:nslot],
                                                    UBbox[:nslot], piiUser)
        fracLB, fracUB = bound_hit_fractions(Phatss[:, :nslot], LBboot[:, :nslot],
                                             UBboot[:, :nslot], 1e-3)
        mcLB, mcUB, mcPII = quantile_mcse(Phatss[:, :nslot], 500)
        MCEs = np.std(Phatss[:, :nslot], axis=0, ddof=1) / np.sqrt(M) if M > 1 else np.full(nslot, np.nan)
        if params.composite is not None:
            cb = params.composite(Phatss)
            comp = np.array([float(np.atleast_1d(params.composite(P_hat[None, :]))[0]),
                             mquantile(cb, 0.025), mquantile(cb, 0.975)])
        else:
            comp = np.full(3, np.nan)
        windows.append(WindowResult(
            tstart=i, P=P_hat, fval=fval, Phatss=Phatss, fvals=fvals, param_estims=param_estims,
            param_bootmedian=bootmed, composite=comp, AICc=(AICc, part1, part2, kpar), noise=noise,
            fit_diagnostics=fit.diagnostics, series=series, Ys=Ys,
            pii={"SCI": PIIv, "PII_rangeFrac": rangeFrac, "PII_scale": scale, "PII_box": PIIbox,
                 "PII_ratio": ratio, "fracLB": fracLB, "fracUB": fracUB, "MCSE_LB": mcLB,
                 "MCSE_UB": mcUB, "MCSE_PII": mcPII, "MCEs": MCEs, "LBbox": LBbox,
                 "UBbox": UBbox, "estimOutsideCI": outside}))
        if verbose:
            lab = ", ".join(f"{_safe_label(params.label[k])}={P_hat[k]:.4g} "
                            f"[{q025[k]:.4g}, {q975[k]:.4g}]"
                            for k in range(params.num) if params.fixed[k] == 0)
            print(f"[{caller}] window {i}: {lab}; AICc = {AICc:.2f}")

    # ---------------- tables over windows ------------------------------------
    labels_triplets, colnames = [], []
    for k in range(params.num):
        lab = params.label[k]
        colnames += [lab, f"{lab}_95%CI LB", f"{lab}_95%CI UB"]
        labels_triplets.append(lab)
    for k in vars.fit_index:
        lab = f"X0({k + 1})"
        colnames += [lab, f"{lab}95%CI LB", f"{lab}95%CI UB"]
        labels_triplets.append(lab)
    extra_slots = []
    if method1 in (3, 4, 5):
        colnames += ["alpha", "alpha_95%CI LB", "alpha_95%CI UB"]
        labels_triplets.append("alpha")
        extra_slots.append(params.num + nfit)
    if method1 == 5:
        colnames += ["d", "d_95%CI LB", "d_95%CI UB"]
        labels_triplets.append("d")
        extra_slots.append(params.num + nfit + 1)
    use_rows = list(range(params.num + nfit)) + extra_slots
    times = [w.tstart for w in windows]
    parrows = [[w.tstart] + list(w.param_estims[use_rows].ravel()) for w in windows]
    partab = pd.DataFrame(parrows, columns=["time"] + colnames)
    comptab = pd.DataFrame([[w.tstart, *w.composite,
                             np.log10(w.composite[2] / w.composite[1]) if w.composite[1] > 0 else np.nan]
                            for w in windows],
                           columns=["time", "composite best fit", "composite 95% CI LB",
                                    "composite 95% CI UB", "composite SCI"])
    aictab = pd.DataFrame(AICcs, columns=["time", "AICc", "AICc part1", "AICc part2", "numparams"])
    sfx = (f"-model_name-{model.name}-fixI0-{params.fixI0}-method-{method1}-dist-{dist1}-tstart-"
           f"{tstart1}-tend-{tend1}-calibrationperiod-{windowsize1}-horizon-{forecastingperiod}-"
           f"{caddisease}-{datatype}.csv")
    if save:
        partab.to_csv(os.path.join(outdir, "parameters-rollingwindow" + sfx), index=False)
        comptab.to_csv(os.path.join(outdir, "parameters-composite" + sfx), index=False)
        aictab.to_csv(os.path.join(outdir, "AICc" + sfx), index=False)
        labs = labels_triplets[:nslot]
        for key, prefix in (("SCI", "SCIs-rollingwindow"), ("PII_box", "PIIbox-rollingwindow"),
                            ("PII_ratio", "PIIratio-rollingwindow"),
                            ("PII_rangeFrac", "PIIrangeFrac-rollingwindow"),
                            ("PII_scale", "PIIscale-rollingwindow"),
                            ("fracLB", "PIIfracLB-rollingwindow"), ("fracUB", "PIIfracUB-rollingwindow")):
            pd.DataFrame([[w.tstart, *w.pii[key]] for w in windows],
                         columns=["time"] + labs).to_csv(os.path.join(outdir, prefix + sfx), index=False)
        mc_cols, mc_rows = ["time"], []
        for lab in labs:
            mc_cols += [f"{lab}_MCSE_CI_LB", f"{lab}_MCSE_CI_UB", f"{lab}_MCSE_PII", f"{lab}_SE_bootmean"]
        mc_cols += ["M_bootstrap", "B2_resamples"]
        for w in windows:
            row = [w.tstart]
            for k in range(nslot):
                row += [w.pii["MCSE_LB"][k], w.pii["MCSE_UB"][k], w.pii["MCSE_PII"][k], w.pii["MCEs"][k]]
            mc_rows.append(row + [M, 500])
        pd.DataFrame(mc_rows, columns=mc_cols).to_csv(os.path.join(outdir, "MCSEs-rollingwindow" + sfx),
                                                      index=False)
        if vars.num > 1:
            last = windows[-1]
            cols = {"time": timevect2}
            for s in range(vars.num):
                Y = last.Ys[s]
                cols[f"median {vars.label[s]}"] = mquantile(Y, 0.5, axis=1)
                cols[f"95% LB {vars.label[s]}"] = mquantile(Y, 0.025, axis=1)
                cols[f"95% UB {vars.label[s]}"] = mquantile(Y, 0.975, axis=1)
            pd.DataFrame(cols).to_csv(os.path.join(outdir, "StateVars" + sfx), index=False)
        np.savez_compressed(
            os.path.join(outdir, f"parameters-ODEModel-{cad}-model_name-{model.name}-fixI0-"
                         f"{params.fixI0}-method-{method1}-dist-{dist1}-tstart-{tstart1}-tend-{tend1}"
                         f"-calibrationperiod-{windowsize1}-forecastingperiod-{forecastingperiod}.npz"),
            param_estims=np.stack([w.param_estims for w in windows], axis=2),
            param_bootmedian=np.column_stack([w.param_bootmedian for w in windows]),
            estimOutsideCI=np.vstack([w.pii["estimOutsideCI"] for w in windows]),
            **{k: np.vstack([w.pii[k] for w in windows]) for k in
               ("PII_box", "PII_ratio", "PII_rangeFrac", "PII_scale", "fracLB", "fracUB",
                "MCSE_LB", "MCSE_UB", "MCSE_PII")})
        if printscreen1:
            # parameter histograms as CSV, as the MATLAB runner does when printscreen1 = 1
            for w in windows:
                for k in range(params.num):
                    if params.fixed[k] == 0:
                        counts, edges = np.histogram(w.Phatss[:, k], bins=10)
                        centers = (edges[:-1] + edges[1:]) / 2
                        pd.DataFrame({"BinEdges": np.round(centers, 4), "Counts": counts}).to_csv(
                            os.path.join(outdir, f"{_safe_label(params.label[k])}-histogram-rollingwindow-"
                                         f"model_name-{model.name}-fixI0-{params.fixI0}-method-{method1}"
                                         f"-dist-{dist1}-tstart-{w.tstart}-calibrationperiod-{windowsize1}"
                                         f"-horizon-{forecastingperiod}-{caddisease}-{datatype}.csv"),
                            index=False)

    requested = M * len(windows)
    if verbose:
        print(f"\nRequested bootstrap datasets: {requested}\nAccepted bootstrap datasets: {accepted}\n"
              f"Unresolved bootstrap datasets: {requested - accepted}\n")
    result = RunResult(options=opts, tstart1=tstart1, tend1=tend1, windowsize1=windowsize1,
                       forecastingperiod=forecastingperiod, windows=windows, AICcs=aictab,
                       parameters=partab, composite=comptab, output_dir=outdir,
                       accepted_bootstrap=accepted, requested_bootstrap=requested)
    if printscreen1:
        from . import plotting
        if forecastingperiod > 0:
            plotting.plot_forecast(result, save=save)
        else:
            plotting.plot_fit(result, save=save)
    return result


def run_forecasting(opts: Options, tstart1=None, tend1=None, windowsize1=None,
                    forecastingperiod=None, *, seed=None, rng=None, n_jobs=1, verbose=True,
                    plot=None, save=True) -> RunResult:
    """Run_Forecasting_ODEModel(options, tstart1, tend1, windowsize1, forecastingperiod).

    Arguments left as None take the value in ``opts``. ``seed`` replaces
    MATLAB's ``rng(seed, 'twister')`` call before the run; ``n_jobs`` runs the
    bootstrap refits in parallel (joblib). ``plot`` overrides printscreen1.
    """
    return _run_windows(opts,
                        opts.tstart1 if tstart1 is None else tstart1,
                        opts.tend1 if tend1 is None else tend1,
                        opts.windowsize1 if windowsize1 is None else windowsize1,
                        opts.forecastingperiod if forecastingperiod is None else forecastingperiod,
                        bool(opts.getperformance), "Run_Forecasting_ODEModel", seed=seed, rng=rng,
                        n_jobs=n_jobs, verbose=verbose, plot=plot, save=save)


def run_fit(opts: Options, tstart1=None, tend1=None, windowsize1=None, *, seed=None, rng=None,
            n_jobs=1, verbose=True, plot=None, save=True) -> RunResult:
    """Run_Fit_ODEModel(options, tstart1, tend1, windowsize1): fit and bootstrap
    the calibration rows only (forecasting period 0)."""
    return _run_windows(opts,
                        opts.tstart1 if tstart1 is None else tstart1,
                        opts.tend1 if tend1 is None else tend1,
                        opts.windowsize1 if windowsize1 is None else windowsize1,
                        0, True, "Run_Fit_ODEModel", seed=seed, rng=rng, n_jobs=n_jobs,
                        verbose=verbose, plot=plot, save=save)


def load_forecast(path):
    """Load a saved Forecast-ODEModel-*.npz file into a dict."""
    with np.load(path, allow_pickle=False) as z:
        out = {k: z[k] for k in z.files}
    if "meta" in out:
        out["meta"] = json.loads(str(out["meta"]))
    return out
