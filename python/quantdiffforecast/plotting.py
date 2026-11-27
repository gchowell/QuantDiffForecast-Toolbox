"""Figures (matplotlib): the plotting halves of Run_Fit_ODEModel.m,
Run_Forecasting_ODEModel.m, plotFit_ODEModel.m, plotForecast_ODEModel.m,
plotODEModel.m and the profile-likelihood figure.

MATLAB drew to the screen with tiledlayout/nexttile; here each function
builds a matplotlib Figure, saves it as a PNG in the output folder and returns
it. Call ``matplotlib.pyplot.show()`` afterwards to see the figures
interactively, or pass ``show=True``.

In MATLAB the plotters re-read the .mat files written by the runners. Here
they take the RunResult returned by run_fit / run_forecasting, so "run, then
plot" needs no file-name bookkeeping:

    res = run_forecasting(opts, 1, 1, 17, 10, seed=1, plot=False)
    plot_forecast(res)
"""
from __future__ import annotations

import os
import re

import numpy as np

from .scoring import mquantile

# one visual system for every figure
INK = "#191919"          # observations, text
MUTED = "#666666"        # axes text, reference lines
GRID = "#e6e6e6"
BLUE = "#257de0"         # the model: median, prediction intervals
ORANGE = "#f47230"       # second series (baseline, point estimate marker)
AQUA = "#00b894"
BOOT = "#b8b8b8"         # individual bootstrap curves


def _plt():
    import matplotlib.pyplot as plt
    return plt


def _style(ax):
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c8c8c8")
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)


def _tex(label):
    """'\\beta' -> '$\\beta$' so matplotlib renders MATLAB-style TeX labels."""
    s = str(label)
    if "$" in s:
        return s
    if "\\" in s or "_" in s or "^" in s:
        return f"${s}$"
    return s


def _fname(s):
    return re.sub(r'[\\/:*?"<>|\s]+', "_", s)


def _finish(fig, path, save, show):
    fig.tight_layout()
    if save and path:
        fig.savefig(path, dpi=160, facecolor="white")
    if show:
        _plt().show()
    else:
        _plt().close(fig)
    return fig


def _ylabel(result, s):
    v = result.options.vars
    j = list(v.fit_index).index(s.fit_index)
    name = v.label[s.fit_index] + ("'(t)" if v.fit_diff[j] else "(t)")
    return f"{name}  ({result.options.caddisease} {result.options.datatype})"


def _fan(ax, s, horizon):
    """Predictive median, 50% and 95% intervals, bootstrap fits and data."""
    t2 = s.timevect2
    f2 = s.forecast2
    q = np.maximum(mquantile(f2, [0.025, 0.25, 0.5, 0.75, 0.975], axis=1), 0)
    ax.fill_between(t2, q[0], q[4], color=BLUE, alpha=0.16, linewidth=0, label="95% prediction interval")
    ax.fill_between(t2, q[1], q[3], color=BLUE, alpha=0.30, linewidth=0, label="50% prediction interval")
    n_show = min(s.fit1.shape[1], 60)
    ax.plot(s.timevect1, s.fit1[:, :n_show], color=BOOT, linewidth=0.5, alpha=0.6)
    ax.plot([], [], color=BOOT, linewidth=1, label="bootstrap model fits")
    ax.plot(t2, q[2], color=BLUE, linewidth=2, label="predictive median")
    nshow = min(s.datalatest.shape[0], t2.size)
    nC = s.timevect1.size
    ax.plot(s.datalatest[:nC, 0], s.datalatest[:nC, 1], "o", color=INK, markersize=5,
            markeredgecolor="white", markeredgewidth=0.8, label="observed (calibration)", zorder=5)
    if horizon > 0 and nshow > nC:
        ax.plot(s.datalatest[nC:nshow, 0], s.datalatest[nC:nshow, 1], "o", markersize=5,
                markerfacecolor="white", markeredgecolor=INK, markeredgewidth=1.2,
                label="observed (held out)", zorder=5)
    if horizon > 0:
        ax.axvline(s.timevect1[-1], color=MUTED, linestyle="--", linewidth=1)
        ax.annotate("forecast start", (s.timevect1[-1], 1), xycoords=("data", "axes fraction"),
                    xytext=(-4, -4), textcoords="offset points", ha="right", va="top",
                    fontsize=8, color=MUTED)
    ax.set_xlim(t2[0], t2[-1])
    ax.set_ylim(0, max(q[4].max() * 1.15, 1))
    ax.set_xlabel("Time")
    _style(ax)


def plot_forecast(result, save=True, show=False):
    """Fit/forecast fan chart, forecast-performance panels, bootstrap
    parameter histograms and state trajectories. Returns the figures."""
    plt = _plt()
    o = result.options
    figs = []
    H = result.forecastingperiod
    for w in result.windows:
        for s in w.series:
            fig, ax = plt.subplots(figsize=(8, 4.8))
            _fan(ax, s, H)
            ax.set_ylabel(_ylabel(result, s))
            kind = "forecast" if H > 0 else "fit"
            ax.set_title(f"{o.model.name}: {kind}, calibration rows {w.tstart}-"
                         f"{w.tstart + result.windowsize1 - 1}", loc="left", fontsize=11, color=INK)
            ax.legend(frameon=False, fontsize=8, loc="upper left")
            tag = _fname(f"{o.model.name}-fit_index-{s.fit_index + 1}-tstart-{w.tstart}-"
                         f"calibrationperiod-{result.windowsize1}-horizon-{H}")
            figs.append(_finish(fig, os.path.join(result.output_dir, f"fig-{kind}-{tag}.png"),
                                save, show))
            if H > 0 and s.performance_forecasting is not None:
                pf = s.performance_forecasting
                fig, axes = plt.subplots(2, 2, figsize=(8, 5.6), sharex=True)
                for ax, col, lab in zip(axes.ravel(), ["MAE", "MSE", "Coverage 95%PI", "WIS"],
                                        ["MAE", "MSE", "Coverage of the 95% PI (%)", "WIS"]):
                    ax.plot(pf["forecasting_horizon"], pf[col], "-o", color=BLUE, linewidth=2,
                            markersize=5)
                    ax.set_title(lab, loc="left", fontsize=10, color=INK)
                    if col == "Coverage 95%PI":
                        ax.axhline(95, color=MUTED, linestyle="--", linewidth=1)
                        ax.set_ylim(0, 105)
                    _style(ax)
                for ax in axes[1]:
                    ax.set_xlabel("Forecasting horizon (mean over leads 1..h)")
                figs.append(_finish(fig, os.path.join(result.output_dir, f"fig-performance-{tag}.png"),
                                    save, show))
        figs.append(plot_parameters(result, w, save=save, show=show))
        if o.vars.num > 1:
            figs.append(plot_states(result, w, save=save, show=show))
    return [f for f in figs if f is not None]


def plot_fit(result, save=True, show=False):
    """plotFit_ODEModel: same panels for a fit without a forecast."""
    return plot_forecast(result, save=save, show=show)


def plot_parameters(result, window=None, save=True, show=False):
    """Bootstrap distributions of the estimated parameters (and the composite)."""
    plt = _plt()
    o = result.options
    w = window or result.windows[-1]
    p = o.params
    items = [(k, _tex(p.label[k]), w.Phatss[:, k], w.P[k], w.param_estims[k, 1:])
             for k in range(p.num) if p.fixed[k] == 0]
    nfit = len(o.vars.fit_index)
    if o.method1 in (3, 4, 5):
        k = p.num + nfit
        items.append((k, r"$\alpha$ (dispersion)", w.Phatss[:, k], w.P[k], w.param_estims[k, 1:]))
    if o.method1 == 5:
        k = p.num + nfit + 1
        items.append((k, "d (variance exponent)", w.Phatss[:, k], w.P[k], w.param_estims[k, 1:]))
    if p.composite is not None:
        cb = np.asarray(p.composite(w.Phatss), float)
        items.append((-1, _tex(p.composite_name or "composite"), cb, w.composite[0], w.composite[1:]))
    if not items:
        return None
    ncols = min(4, len(items))
    nrows = int(np.ceil(len(items) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 3.0 * nrows), squeeze=False)
    for ax, (_, lab, draws, est, ci) in zip(axes.ravel(), items):
        ax.hist(draws, bins=20, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.8)
        ax.axvline(est, color=ORANGE, linestyle="--", linewidth=2)
        ax.set_xlabel(lab)
        ax.set_title(f"{est:.3g}\n95% CI: {ci[0]:.3g}, {ci[1]:.3g}", fontsize=9, color=INK, loc="left")
        _style(ax)
    for ax in axes.ravel()[len(items):]:
        ax.set_visible(False)
    axes[0, 0].set_ylabel("Bootstrap frequency")
    fig.suptitle(f"{o.model.name}: bootstrap estimates (dashed line = fit to the observed data), "
                 f"window {w.tstart}", fontsize=10, color=INK, x=0.01, ha="left")
    tag = _fname(f"{o.model.name}-tstart-{w.tstart}-calibrationperiod-{result.windowsize1}")
    return _finish(fig, os.path.join(result.output_dir, f"fig-parameters-{tag}.png"), save, show)


def plot_states(result, window=None, save=True, show=False):
    """Bootstrap median and 95% band of every model state."""
    plt = _plt()
    o = result.options
    w = window or result.windows[-1]
    n = o.vars.num
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    t2 = w.series[0].timevect2
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.9 * nrows), squeeze=False, sharex=True)
    for s, ax in zip(range(n), axes.ravel()):
        q = mquantile(w.Ys[s], [0.025, 0.5, 0.975], axis=1)
        ax.fill_between(t2, q[0], q[2], color=BLUE, alpha=0.2, linewidth=0)
        ax.plot(t2, q[1], color=BLUE, linewidth=2)
        ax.set_title(_tex(o.vars.label[s]), loc="left", fontsize=10, color=INK)
        _style(ax)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Time")
    fig.suptitle(f"{o.model.name}: model states, bootstrap median and 95% band (window {w.tstart})",
                 fontsize=10, color=INK, x=0.01, ha="left")
    tag = _fname(f"{o.model.name}-tstart-{w.tstart}-calibrationperiod-{result.windowsize1}")
    return _finish(fig, os.path.join(result.output_dir, f"fig-states-{tag}.png"), save, show)


def plot_ode_model(opts, windowsize1=None, M=20, seed=1, save=True, show=False):
    """plotODEModel: explore model trajectories before fitting.

    Draws M parameter vectors (free parameters uniform within their bounds; the
    first draw is params.initial), integrates the model over the calibration
    window and marks the draw with the smallest sum of squared errors.
    """
    from .config import load_data, make_window, validate_options
    from .core import initial_conditions, observation_curve, solve_ode
    plt = _plt()
    validate_options(opts, "plotODEModel")
    data = load_data(opts)
    W = windowsize1 or opts.windowsize1
    i0 = opts.tstart1 - 1
    d = data[i0:i0 + W]
    p, v = opts.params, opts.vars
    v.window = make_window(data, opts.tstart1, W)
    rng = np.random.default_rng(seed)
    nfit = len(v.fit_index)
    curves, sses, states = [], [], []
    for j in range(M):
        for _ in range(50):
            pj = p.initial.copy()
            if j > 0:
                free = p.fixed == 0
                pj[free] = rng.uniform(p.LB[free], p.UB[free])
            z = np.concatenate([pj, d[0, 1:], [1, 1]])
            try:
                F = solve_ode(opts.model.fc, d[:, 0], initial_conditions(z, p, v), z, p.extra0)
                c = observation_curve(F, v).reshape(W, nfit, order="F")
            except Exception:
                if j == 0:
                    break
                continue
            curves.append(c)
            states.append(F)
            sses.append(float(np.sum((c - d[:, 1:]) ** 2)))
            break
    if not curves:
        raise RuntimeError("No parameter draw could be integrated.")
    best = int(np.argmin(sses))
    fig, axes = plt.subplots(1, nfit, figsize=(6.2 * nfit, 4.4), squeeze=False)
    for j, ax in enumerate(axes[0]):
        for c in curves:
            ax.plot(d[:, 0], c[:, j], color=BOOT, linewidth=0.8)
        ax.plot([], [], color=BOOT, label="sampled parameter values")
        ax.plot(d[:, 0], curves[best][:, j], color=AQUA, linewidth=2.2, label="smallest SSE among the draws")
        ax.plot(d[:, 0], d[:, j + 1], "o", color=INK, markersize=5, markeredgecolor="white",
                label="observed")
        ymax = max(d[:, j + 1].max() * 3, 1)
        ax.set_ylim(0, ymax)
        ax.set_xlabel("Time")
        ax.set_ylabel(v.label[v.fit_index[j]] + ("'(t)" if v.fit_diff[j] else "(t)"))
        ax.legend(frameon=False, fontsize=8)
        _style(ax)
    fig.suptitle(f"{opts.model.name}: trajectories before fitting", fontsize=10, color=INK,
                 x=0.01, ha="left")
    return _finish(fig, os.path.join(opts.outdir(), f"fig-model-{_fname(opts.model.name)}.png"),
                   save, show)


def plot_profiles(PL, output_dir=None, save=True, show=False):
    """Profile-likelihood figure: deviance against each profiled parameter."""
    plt = _plt()
    nk = len(PL["param"])
    ncols = min(nk, 4)
    nrows = int(np.ceil(nk / ncols))
    thr = PL["threshold"]
    level = 100 * (1 - PL["alpha"])
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.8 * nrows), squeeze=False)
    for p, ax in zip(PL["param"], axes.ravel()):
        ok = np.isfinite(p["D"])
        ax.plot(p["theta"][ok], p["D"][ok], "-o", color=BLUE, linewidth=1.8, markersize=4,
                label="profile deviance")
        ax.axhline(thr, color=ORANGE, linestyle="--", linewidth=1.3,
                   label=f"{level:.0f}% threshold")
        ci = p["CI"]
        if np.all(np.isfinite(ci)):
            ax.axvspan(ci[0], ci[1], color=BLUE, alpha=0.10, linewidth=0,
                       label=f"{level:.0f}% interval")
        ax.plot([p["mle"]], [0], "D", color=INK, markersize=6, label="MLE", zorder=5)
        th = p["theta"][ok]
        if th.size and th.min() > 0 and th.max() / th.min() > 10:
            ax.set_xscale("log")
        ax.set_ylim(0, max(1.5 * thr, min(np.nanmax(np.append(p["D"][ok], 0)) * 1.1, 6 * thr)))
        ax.set_xlabel(_tex(p["label"]))
        fmt = lambda x: "NaN" if not np.isfinite(x) else f"{x:.4g}"
        ax.set_title(f"MLE {p['mle']:.4g};  {level:.0f}% PL interval [{fmt(ci[0])}, {fmt(ci[1])}]\n"
                     f"{p['status']}", fontsize=8.5, color=INK, loc="left")
        _style(ax)
    for ax in axes.ravel()[nk:]:
        ax.set_visible(False)
    axes[0, 0].set_ylabel(r"$-2\,\Delta \log L$")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{PL['model']}: profile likelihood, rows {PL['tstart1']}-"
                 f"{PL['tstart1'] + PL['windowsize1'] - 1}, method {PL['method1']}",
                 fontsize=10, color=INK, x=0.01, ha="left")
    path = None
    if output_dir:
        path = os.path.join(output_dir, "PLprofiles-" + PL.get("suffix", "profile").replace(".csv", "") + ".png")
    return _finish(fig, path, save, show)


def plot_ensemble(members, labels, weights, ens_summary, ens_label, timevect1, timevect2,
                  datalatest, horizon, path=None, save=True, show=False):
    """One panel per member plus the ensemble (run_ensemble's figure)."""
    plt = _plt()
    K = len(members)
    fig, axes = plt.subplots(1, K + 1, figsize=(4.6 * (K + 1), 4.2), squeeze=False, sharey=True)
    panels = [(labels[k] + f" (w = {weights[k]:.2f})",
               mquantile(members[k], [0.5, 0.025, 0.975], axis=1).T) for k in range(K)]
    panels.append((ens_label, ens_summary))
    nC = len(timevect1)
    for ax, (ttl, q) in zip(axes[0], panels):
        ax.fill_between(timevect2, np.maximum(q[:, 1], 0), q[:, 2], color=BLUE, alpha=0.18, linewidth=0,
                        label="95% prediction interval")
        ax.plot(timevect2, q[:, 0], color=BLUE, linewidth=2, label="median")
        n = min(datalatest.shape[0], len(timevect2))
        ax.plot(datalatest[:nC, 0], datalatest[:nC, 1], "o", color=INK, markersize=4.5,
                markeredgecolor="white", label="observed (calibration)")
        if n > nC:
            ax.plot(datalatest[nC:n, 0], datalatest[nC:n, 1], "o", markersize=4.5,
                    markerfacecolor="white", markeredgecolor=INK, label="observed (held out)")
        if horizon > 0:
            ax.axvline(timevect1[-1], color=MUTED, linestyle="--", linewidth=1)
        ax.set_title(ttl, loc="left", fontsize=10, color=INK)
        ax.set_xlabel("Time")
        _style(ax)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="upper left")
    return _finish(fig, path, save, show)
